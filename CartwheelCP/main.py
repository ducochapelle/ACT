# Needs. More. OOD.

# add option to add vertices. 
#   option: You can even do it cool where if centernode is edge node, pick one from vertex pool.
# add warning if midnode is in edgenodes
# add on activate scope put selection in selection manager
# make default spoke diameter yellow
# make spoke diameter automatic
# add sideways force overdracht
#   shaft hole beam link link triangle?
#   hole hole
# Evade rotating nodes of model by either using 'CE in direction' or 'double spokes (----||****O****||----) where ( = edge, - = link, || = connection in plane, * = beam, 0 = shaft'
# add unit of density, diameter and emodulus properly (don't guess, ask ansys!)
# D uit CP in voor constrain in Z HIGH PRIORITY


from datetime import datetime
drawingObjects = {}
NearestNodeIdByEdgeId = {}

def ExtAPILogWriteMessage(string):
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+string)
    pass
    
def ButtonClick1(analysis):
    load = analysis.CreateLoadObject("CartwheelCP")
                
def ShowCartwheel(load):
    def getRealCenterNode(nodeXYZ,load):
        MeshData = load.Analysis.MeshData
        Vertex = load.Properties.GetByName("Vertex").Value
        LocationByNodeId = {}
        while True:
            if Vertex in ("Select Vertex"):
                VertexIds = load.Properties.GetByName("Vertices").Value.Ids
                NodeIds = [MeshData.MeshRegionById(VertexId).NodeIds[0] for VertexId in VertexIds]
            elif Vertex in ("Automatic Vertex"):
                NodeIds = MeshData.NodeIds
            elif Vertex in ("Create Vertex"):
                break
            for NodeId in NodeIds:
                LocationByNodeId[NodeId] = ( MeshData.NodeById(NodeId).X,MeshData.NodeById(NodeId).Y,MeshData.NodeById(NodeId).Z )
            break
        ExtAPILogWriteMessage(str(Vertex)+" "+str(LocationByNodeId))    
        distances = {}
        for node in LocationByNodeId:
            distances[node] = ( (LocationByNodeId[node][0] - nodeXYZ[0])**2 + (LocationByNodeId[node][1] - nodeXYZ[1])**2 + (LocationByNodeId[node][2] - nodeXYZ[2])**2 )**(0.5)
        return min(distances, key=distances.get)

    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+"START Cartwheel display")

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
    
    ExtAPILogWriteMessage("Geom is: "+str(Edges.IsValid)+"...")

    rainbow = {0:0xFF0000
              ,1:0xFF6600
              ,2:0xFFFF00
              ,3:0x99FF00
              ,4:0x00FF00
              ,5:0x00CCFF
              ,6:0x0099FF
              ,7:0x3333FF
              ,8:0x990099} ; i = 0
    
    if Edges.IsValid and MeshData.ElementCount:
        ExtAPILogWriteMessage("1...")
        drawingObjects[load.Id] = ExtAPI.Graphics.CreateAndOpenDraw3DContext()
        drawingObjects[load.Id].Color = load.Color
        drawingObjects[load.Id].LineWeight = 2
        pts = System.Array.CreateInstance(float,6)   
        EdgeIds = Edges.Value.Ids
        for EdgeId in EdgeIds:
            ExtAPILogWriteMessage("2...")
            NodeIds =   MeshData.MeshRegionById(EdgeId).NodeIds
            NodeCount = MeshData.MeshRegionById(EdgeId).NodeCount
            # virtual centernode
            pts[0] = sum([MeshData.NodeById(n).X for n in NodeIds],0.0)/NodeCount
            pts[1] = sum([MeshData.NodeById(n).Y for n in NodeIds],0.0)/NodeCount
            pts[2] = sum([MeshData.NodeById(n).Z for n in NodeIds],0.0)/NodeCount
            if Vertex in ("Automatic Vertex", "Select Vertex"):
                # real centernode
                if not EdgeId in NearestNodeIdByEdgeId:
                    NearestNodeIdByEdgeId[EdgeId] = getRealCenterNode(pts,load)
                nearestNode = NearestNodeIdByEdgeId[EdgeId]
                ExtAPILogWriteMessage("3... "+str(nearestNode))
                pts[0] = MeshData.NodeById(nearestNode).X
                pts[1] = MeshData.NodeById(nearestNode).Y
                pts[2] = MeshData.NodeById(nearestNode).Z
            for NodeId in MeshData.MeshRegionById(EdgeId).NodeIds:
                ExtAPILogWriteMessage("4...")
                pts[3] = MeshData.NodeById(NodeId).X
                pts[4] = MeshData.NodeById(NodeId).Y
                pts[5] = MeshData.NodeById(NodeId).Z
                if load.Properties.GetByName("Colors").Value == "Rainbow":
                    drawingObjects[load.Id].Color = rainbow[i] ; i += 1 ; i = i % rainbow.Count
                drawingObjects[load.Id].DrawPolyline(pts)
        drawingObjects[load.Id].Close()
        drawingObjects[load.Id].Visible = True
    ExtAPILogWriteMessage(str(drawingObjects))
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+"END Cartwheel display")

def HideCartwheel(load):
    global drawingObjects
    if load.Id in drawingObjects:
        if drawingObjects[load.Id]!=None:
            drawingObjects[load.Id].Visible = False
        drawingObjects[load.Id] = None
    pass
    
def CartwheelCP(load, stream):    
    ExtAPILogWriteMessage("Cartwheel... 666!")
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
    stream.Write("C*** halfway cartwheels \n")
    stream.Write("TYPE,zSpokeEtype \n")
    stream.Write("REAL,zSpokeEtype \n")


    for geomId in geomIds:
        meshRegion = mesh.MeshRegionById(geomId)
        nodeIds = meshRegion.NodeIds
        
        #align coordinate system on first three edge nodes
        stream.Write("/PREP7 \n")
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
                stream.Write("NSEL,a,,,shaftNode"+geomId.ToString()+" \n")
                stream.Write("CP,NEXT,UZ,shaftNode"+geomId.ToString()+",lastNode \n")
                stream.Write("*GET,zNset,CP,0,MAX \n")
                token = True
            else:
                stream.Write("C***, handenopelkaar \n")
                stream.Write("CP,zNset,,lastNode \n")
        token = False
        
        #additional displacements
        stream.Write("NSEL,,,,shaftNode"+geomId.ToString()+" \n")
        stream.Write("CSYS,CSmax+1 \n")        
        stream.Write("NROTAT,all \n")
        stream.Write("CSYS,0 \n")
        stream.Write("/SOLU \n")
        if Vertex == "Create Vertex":
            stream.Write("D,ALL,ALL \n")
        else:
            if load.Properties.GetByName("ConstrainShaftU").Value == "Yes":
                stream.Write("D,ALL,UZ,0 \n")   # beam line = Z <--------------------------------------------------------
            if load.Properties.GetByName("ConstrainShaftRot").Value == "Yes":
                stream.Write("D,ALL,ROTZ,0 \n") # beam line = Z <--------------------------------------------------------
        stream.Write("/PREP7 \n")
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
        load.Properties.GetByName("ConstrainShaftU").Value = "No"
        load.Properties.GetByName("ConstrainShaftRot").Value = "No"
    elif prop.Value == "Select Vertex":
        load.Properties.GetByName("Vertices").Visible = True
        load.Properties.GetByName("ConstrainShaftU").ReadOnly = False
        load.Properties.GetByName("ConstrainShaftRot").ReadOnly = False
        load.Properties.GetByName("ConstrainShaftU").Value = "No"
        load.Properties.GetByName("ConstrainShaftRot").Value = "No"
    elif prop.Value == "Create Vertex":
        load.Properties.GetByName("Vertices").Visible = False
        load.Properties.GetByName("ConstrainShaftU").ReadOnly = True
        load.Properties.GetByName("ConstrainShaftRot").ReadOnly = True
        load.Properties.GetByName("ConstrainShaftU").Value = "Yes"
        load.Properties.GetByName("ConstrainShaftRot").Value = "Yes"
    HideCartwheel(load)
    ShowCartwheel(load)

def SelectColors(load,prop):
    prop.ClearOptions()
    prop.AddOption("Default")
    prop.AddOption("Rainbow")
    HideCartwheel(load)
    ShowCartwheel(load)
    
def SelectNopr(load,prop):
    prop.ClearOptions()
    prop.AddOption("/NOPR")
    prop.AddOption("/GOPR")
    
