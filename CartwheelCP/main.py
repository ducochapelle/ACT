# Needs. More. OOD.
# Maybe the bug out is when mechanical is not closed for 48+ hours?
# Maybe the bug out can be prevented by FIRST checking the ACT plugins in the project manager and THEN opening the project

# make refresh mesh work when you clicked on this stuff when you didnt even had mesh
# add color option RGB
# add warning if midnode is in edgenodes
# add on activate scope put selection in selection manager
# make default spoke diameter yellow
# make spoke diameter automatic
# add sideways force overdracht
#   shaft hole beam link link triangle?
#   hole hole
# Evade rotating nodes of model by either using 'CE in direction' or 'double spokes (----||****O****||----) where ( = edge, - = link, || = connection in plane, * = beam, 0 = shaft'
# add unit of density, diameter and emodulus properly (don't guess, ask ansys!)
# getRealCenterNode moet sneller
# Less flood with /nopr
# when suppressed -> not visible
# ask area instead of diameter (or maybe even auto-calculate with edge-area-thickness, line length and node count)
# change spoke diameter to spoke area
if ExtAPI.Context == "Mechanical":
    from math import hypot
    from datetime import datetime
    drawingObjects = {}
    NearestNodeIdByEdgeId = {}
    GlobalLocationByNodeId = {}
    global ColorTable
    ColorTable =    {"Red"       :0xFF0000,
                     "Maroon"    :0x800000,
                     "Brown"     :0xA52A2A,
                     "Orange"    :0xFFA500,
                     "Yellow"    :0xFFFF00,
                     "Lime"      :0x00FF00,
                     "Green"     :0x008000,
                     "Olive"     :0x808000,
                     "Cyan"      :0x00FFFF,
                     "Blue"      :0x0000FF,
                     "LightBlue" :0xADD8E6,
                     "DarkBlue"  :0x0000A0,
                     "Fuchsia"   :0xFF00FF,
                     "Purple"    :0x800080,
                     "White"     :0xFFFFFF,
                     "Silver"    :0xC0C0C0,
                     "Grey"      :0x808080,
                     "Black"     :0x000000}
    
def ExtAPILogWriteMessage(string):
    # ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+string)
    pass
    
def ButtonClick1(analysis):
    load = analysis.CreateLoadObject("CartwheelCP")
                
def getRealCenterNode(location_in_space,load, Vertex, MeshData):
    ExtAPILogWriteMessage("getRealCenterNode 1")
    LocationByNodeId = {}
    if Vertex == "Select Vertex":
        VertexIds = load.Properties.GetByName("Vertices").Value.Ids
        NodeIds = [MeshData.MeshRegionById(VertexId).NodeIds[0] for VertexId in VertexIds]
        for NodeId in NodeIds:
            node = MeshData.NodeById(NodeId)
            LocationByNodeId[NodeId] = (node.X,node.Y,node.Z)
    elif Vertex == "Automatic Vertex":
        global GlobalLocationByNodeId
        if GlobalLocationByNodeId == {}:
            NodeIds = MeshData.NodeIds
            for NodeId in NodeIds:
                node = MeshData.NodeById(NodeId)
                GlobalLocationByNodeId[NodeId] = (node.X,node.Y,node.Z)
        LocationByNodeId = GlobalLocationByNodeId
        
    # http://stackoverflow.com/questions/16979618/given-a-set-of-locations-and-a-single-location-find-the-closest-location-from-t
    # http://en.wikipedia.org/wiki/Nearest_neighbor_search
    # http://en.wikipedia.org/wiki/R-tree
    closest_node = None
    closest_distance = 1e100  # An arbitrary, HUGE, value
    x,y,z = location_in_space[:3]
    for NodeId, NodeLocation in LocationByNodeId.iteritems():
        distance = (NodeLocation[0] - x)**2 + (NodeLocation[1] - y)**2 + (NodeLocation[2] - z)**2
        if distance < closest_distance:
            closest_distance = distance
            closest_node = NodeId
    return closest_node


def ShowCartwheel(load):
    if load.Analysis.MeshData.NodeCount == 0:
        return
    global NearestNodeIdByEdgeId
    global drawingObjects
    if not load.Id in drawingObjects:
        drawingObjects[load.Id] = None
    if drawingObjects[load.Id] != None:
        drawingObjects[load.Id].Visible = False     
    Edges = load.Properties.GetByName("One")
    Vertices = load.Properties.GetByName("Vertices")
    Vertex = load.Properties.GetByName("Vertex").Value
    if Vertex == "Select Vertex" and not Vertices.IsValid:
        return
        
    MeshData = load.Analysis.MeshData

    rainbow = {0:0xFF0000
              ,1:0xFF6600
              ,2:0xFFFF00
              ,3:0x99FF00
              ,4:0x00FF00
              ,5:0x00CCFF
              ,6:0x0099FF
              ,7:0x3333FF
              ,8:0x990099} ; i = 0
    global ColorTable
    
    if Edges.IsValid and MeshData.ElementCount:
        drawingObjects[load.Id] = ExtAPI.Graphics.CreateAndOpenDraw3DContext()
        drawingObjects[load.Id].Color = load.Color
        drawingObjects[load.Id].LineWeight = 2
        Color = load.Properties.GetByName("Colors").Value
        pts = System.Array.CreateInstance(float,6)   
        EdgeIds = Edges.Value.Ids
        for EdgeId in EdgeIds:
            NodeIds =   MeshData.MeshRegionById(EdgeId).NodeIds
            NodeCount = MeshData.MeshRegionById(EdgeId).NodeCount
            # virtual centernode
            pts[0] = sum([MeshData.NodeById(n).X for n in NodeIds],0.0)/NodeCount
            pts[1] = sum([MeshData.NodeById(n).Y for n in NodeIds],0.0)/NodeCount
            pts[2] = sum([MeshData.NodeById(n).Z for n in NodeIds],0.0)/NodeCount
            if Vertex in ("Automatic Vertex", "Select Vertex"):
                # real centernode
                if not EdgeId in NearestNodeIdByEdgeId:
                    NearestNodeIdByEdgeId[EdgeId] = getRealCenterNode(pts,load, Vertex, MeshData)
                    ExtAPILogWriteMessage("getRealCenterNode 3")
                nearestNode = NearestNodeIdByEdgeId[EdgeId]
                pts[0] = MeshData.NodeById(nearestNode).X
                pts[1] = MeshData.NodeById(nearestNode).Y
                pts[2] = MeshData.NodeById(nearestNode).Z
            for NodeId in MeshData.MeshRegionById(EdgeId).NodeIds:
                pts[3] = MeshData.NodeById(NodeId).X
                pts[4] = MeshData.NodeById(NodeId).Y
                pts[5] = MeshData.NodeById(NodeId).Z
                if Color == "Default":
                    pass
                elif Color == "Rainbow":
                    drawingObjects[load.Id].Color = rainbow[i] ; i += 1 ; i = i % rainbow.Count
                elif Color == "Custom":
                    drawingObjects[load.Id].Color = int(load.Properties.GetByName("RGB").Value,16)
                else:
                    drawingObjects[load.Id].Color = ColorTable[Color]
                drawingObjects[load.Id].DrawPolyline(pts)
        drawingObjects[load.Id].Close()
        drawingObjects[load.Id].Visible = True

def HideCartwheel(load):
    global drawingObjects
    if load.Id in drawingObjects:
        if drawingObjects[load.Id]!=None:
            drawingObjects[load.Id].Visible = False
        drawingObjects[load.Id] = None
    pass
    
def CartwheelCP(load, stream):    
    # Get the scoped geometry:
    propGeo = load.Properties.GetByName("One")
    Vertex = load.Properties.GetByName("Vertex").Value
    MeshData = load.Analysis.MeshData
    geomIds = propGeo.Value.Ids
    
    mesh = ExtAPI.DataModel.MeshDataByName("Global")
    
    stream.Write("C*** Applying cartwheels \n")
    stream.Write(load.Properties.GetByName("nopr").Value + "\n")
    stream.Write("/PREP7 \n")
    stream.Write("sstif,on !activate non-lin procedure needed to calculate tension only \n")
    
    #Material properties
    stream.Write("*GET,zMat,MAT,0,NUM,MAX \n")
    stream.Write("zMat = zMat + 1 \n")
    stream.Write("MPTEMP,,                 \n")
    stream.Write("MPTEMP,1,0               \n")
    stream.Write("MPDATA,DENS,zMat,,"+load.Properties.GetByName("Density").Value+" \n")
    stream.Write("MPDATA,EX,zMat,,"+load.Properties.GetByName("Emodulus").Value+" \n")
    stream.Write("MPDATA,PRXY,zMat,,0.3        \n")
    stream.Write("MAT,zMat \n")
    
    #Spoke properties
    spokeDiameter = load.Properties.GetByName("SpokeDiameter")
    stream.Write("ET,,LINK180 \n")
    stream.Write("*GET,zSpokeEtype,ETYP,0,NUM,MAX \n")
    stream.Write("R,zSpokeEtype,"+spokeDiameter.Value+"**2,,-1 \n")
    
    #create spokes
    stream.Write("TYPE,zSpokeEtype \n")
    stream.Write("REAL,zSpokeEtype \n")


    for geomId in geomIds:
        meshRegion = mesh.MeshRegionById(geomId)
        nodeIds = meshRegion.NodeIds
        
        #align coordinate system on first three edge nodes
        # stream.Write("/PREP7 \n")
        stream.Write("CSYS,0 \n")
        stream.Write("NSEL,ALL \n")
        stream.Write("NWPLAN,,"+nodeIds[0].ToString()+","+nodeIds[1].ToString()+","+nodeIds[2].ToString()+" \n")
        
        #offset coordinate system to get center node
        stream.Write("NSEL,NONE \n")
        for nodeId in nodeIds:
            if not nodeId % (len(nodeIds)/9+1): #NWPAVE has limit of 9:
                stream.Write("NSEL,A,,,"+nodeId.ToString()+" \n" )
        stream.Write("NWPAVE, ALL \n")
        stream.Write("CSYS, 4 \n")
        if Vertex == "Automatic Vertex":
            stream.Write("ALLSEL \n")        
            stream.Write("shaftNode"+geomId.ToString()+" = NODE(0,0,0) \n")
        elif Vertex == "Select Vertex":
            stream.Write("NSEL,None \n")
            NodeIds = [MeshData.MeshRegionById(VertexId).NodeIds[0] for VertexId in load.Properties.GetByName("Vertices").Value.Ids]
            for NodeId in NodeIds:
                stream.Write("NSEL,A,,,"+str(NodeId)+"  \n")
            stream.Write("shaftNode"+geomId.ToString()+" = NODE(0,0,0) \n")
        elif Vertex == "Create Vertex":
            stream.Write("NSEL,None \n")
            stream.Write("N \n")
            stream.Write("shaftNode"+geomId.ToString()+" = NDNEXT(0) \n")
            stream.Write("Nlist \n")
            stream.Write("ALLSEL \n")        
            
        stream.Write("*get,CSmax,cdsy,,NUM,MAX $ CSmax = max(CSmax,11) \n")
        stream.Write("CSWPLA,CSmax+1,0 \n")
        
        #create spokes between center node and coupled edge nodes
        token = False
        for nodeId in nodeIds:    
            stream.Write("NSEL,NONE \n")
            stream.Write("NWPAVE,"+nodeId.ToString()+" \n")
            stream.Write("CSYS, 4 \n")
            stream.Write("N \n")
            stream.Write("CSYS, CSmax+1 \n")
            stream.Write("NROTAT,all \n")
            stream.Write("CSYS, 4 \n")
            stream.Write("lastNode = ndnext(0) \n")
            stream.Write("E,lastNode,shaftNode"+geomId.ToString()+" \n")
            stream.Write("NSEL,A,,,"+nodeId.ToString()+" \n")
            stream.Write("CSYS,CSmax+1 \n")
            stream.Write("NROTAT,all \n")
            stream.Write("CP,NEXT,UX,ALL \n")      # beam line = Z <--------------------------------------------------------
            stream.Write("CP,NEXT,UY,ALL \n")       # beam line = Z <--------------------------------------------------------
            stream.Write("CSYS, 4 \n")
            if token == False:
                stream.Write("NSEL,a,,,shaftNode"+geomId.ToString()+" \n") # SHOULD NOT COUPLE nodeId!!! maar het lijkt gewoon te kloppen
                stream.Write("CP,NEXT,UZ,shaftNode"+geomId.ToString()+",lastNode \n") # Kijk dan, hij coupled aleen shaftnode en lastnode, nooit nodeId!
                stream.Write("*GET,zNset,CP,0,MAX \n")
                token = True
            else:
                stream.Write("CP,zNset,,lastNode \n")
        token = False
        
        #additional displacements
        stream.Write("NSEL,,,,shaftNode"+geomId.ToString()+" \n")
        stream.Write("CSYS,CSmax+1 \n")        
        stream.Write("NROTAT,all \n")
        stream.Write("CSYS,0 \n")
        stream.Write("/SOLU \n")
        stream.Write(load.Properties.GetByName("nopr").Value + "\n")    #redundant?
        if Vertex == "Create Vertex":
            stream.Write("D,ALL,ALL \n")
        else:
            if load.Properties.GetByName("ConstrainShaftU").Value == "Yes":
                stream.Write("D,ALL,UZ,0 \n")   # beam line = Z <--------------------------------------------------------
            if load.Properties.GetByName("ConstrainShaftRot").Value == "Yes":
                stream.Write("D,ALL,ROTZ,0 \n") # beam line = Z <--------------------------------------------------------
        stream.Write("/PREP7 \n")
        stream.Write(load.Properties.GetByName("nopr").Value + "\n")    #redundant?
        stream.Write("CSYS,4 \n")
        
    
    stream.Write("ALLS \n")
    stream.Write("CSYS,0 \n")
    stream.Write("zMat = \n")
    stream.Write("zSpokeEtype = \n")
    stream.Write("shaftNode = \n")
    stream.Write("CSmax = \n")
    stream.Write("lastNode = \n")
    stream.Write("/GOPR \n")
    stream.Write("/SOLU \n")
    stream.Write("C*** Finished applying cartwheels \n")

def SelectConstrainShaft(load,prop):
    prop.ClearOptions()
    prop.AddOption("Yes")
    prop.AddOption("No")
    
def SelectVertex(load,prop):
    prop.ClearOptions()
    prop.AddOption("Automatic Vertex")
    prop.AddOption("Select Vertex")
    prop.AddOption("Create Vertex")

def ValidateVertex(load, prop):
    global NearestNodeIdByEdgeId
    NearestNodeIdByEdgeId = {}
    if prop.Value == "Automatic Vertex":
        load.Properties.GetByName("Vertices").Visible = False
        load.Properties.GetByName("ConstrainShaftU").ReadOnly = False
        load.Properties.GetByName("ConstrainShaftRot").ReadOnly = False
        # load.Properties.GetByName("ConstrainShaftU").Value = "No"
        # load.Properties.GetByName("ConstrainShaftRot").Value = "No"
    elif prop.Value == "Select Vertex":
        load.Properties.GetByName("Vertices").Visible = True
        load.Properties.GetByName("ConstrainShaftU").ReadOnly = False
        load.Properties.GetByName("ConstrainShaftRot").ReadOnly = False
        # load.Properties.GetByName("ConstrainShaftU").Value = "No"
        # load.Properties.GetByName("ConstrainShaftRot").Value = "No"
    elif prop.Value == "Create Vertex":
        load.Properties.GetByName("Vertices").Visible = False
        load.Properties.GetByName("ConstrainShaftU").ReadOnly = True
        load.Properties.GetByName("ConstrainShaftRot").ReadOnly = True
        load.Properties.GetByName("ConstrainShaftU").Value = "Yes"
        load.Properties.GetByName("ConstrainShaftRot").Value = "Yes"
    # HideCartwheel(load)
    # ShowCartwheel(load)
    pass

def ClearedMesh(load,prop):
    global GlobalLocationByNodeId
    GlobalLocationByNodeId = {}
    global NearestNodeIdByEdgeId
    NearestNodeIdByEdgeId = {}
    HideCartwheel(load)
    ShowCartwheel(load)

def IsValid(foo,bar):
    return True
    
def SelectColors(load,prop):
    prop.ClearOptions()
    prop.AddOption("Default")
    prop.AddOption("Rainbow")
    prop.AddOption("Custom")
    global ColorTable
    for c in ColorTable:
        prop.AddOption(c)
    HideCartwheel(load)
    ShowCartwheel(load)

def ValidateColors(load,prop):
    load.Properties.GetByName("RGB").Visible = True if prop.Value == "Custom" else False

    
def SelectNopr(load,prop):
    prop.ClearOptions()
    prop.AddOption("/NOPR")
    prop.AddOption("/GOPR")
    
