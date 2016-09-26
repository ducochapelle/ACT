# named seelction selection method does not work

def apdl(string):
    stream.Write(string+ " \n")

def init(context):    
    ExtAPI.Log.WriteMessage("Init APDL...")

def ButtonClick1(analysis):
    ExtAPI.Log.WriteMessage("ButtonClick1...")
    # clr.AddReference("System.Diagnostics.Debugger")
    # import System
    # System.Diagnostics.Launch()
    load = analysis.CreateLoadObject("CartwheelCP")

def CartwheelCP(load, stream):    
    ExtAPI.Log.WriteMessage("Cartwheel... 666!")
    global g_load
    g_load = load
    # Get the scoped geometry:
    propGeo = load.PropertyByName("One")
    geomIds = propGeo.Value
    
    mesh = ExtAPI.DataModel.MeshDataByName("Global")
    
    stream.Write("C*** Applying cartwheels \n")
    stream.Write("/PREP7 \n")
    stream.Write("sstif,on !activate non-lin procedure needed to calculate tension only \n")
    
    #Material properties
    stream.Write("*GET,zMat,MAT,0,NUM,MAX \n")
    stream.Write("zMat = zMat + 1 \n")
    stream.Write("MPTEMP,,                 \n")
    stream.Write("MPTEMP,1,0               \n")
    stream.Write("MPDATA,DENS,zMat,,7.86e-6   \n")
    stream.Write("MPDATA,EX,zMat,,2.1e5       \n")
    stream.Write("MPDATA,PRXY,zMat,,0.3        \n")
    stream.Write("MAT,zMat \n")
    
    #Spoke properties
    spokeDiameter = load.PropertyByName("SpokeDiameter")
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
        stream.Write("NSEL,ALL \n")
        stream.Write("NWPLAN,,"+nodeIds[0].ToString()+","+nodeIds[1].ToString()+","+nodeIds[2].ToString()+" \n")
        
        #offset coordinate system to get center node
        stream.Write("NSEL,NONE \n")
        for nodeId in nodeIds:
            stream.Write("NSEL,A,,,"+nodeId.ToString()+" \n" )
        stream.Write("NWPAVE, ALL \n")
        stream.Write("ALLSEL \n")        
        stream.Write("CSYS, 4 \n")
        stream.Write("shaftNode"+geomId.ToString()+" = NODE(0,0,0) \n")
        stream.Write("*get,CSmax,cdsy,,NUM,MAX $ CSmax = max(CSmax,11) \n")
        stream.Write("CSWPLA,CSmax+1 \n")
        # stream.Write("CLOCAL,CSmax+1,0,,,,,,-90 \n")
        
        #create spokes between center node and coupled edge nodes
        for nodeId in nodeIds:    
            stream.Write("NSEL,NONE \n")
            stream.Write("NWPAVE,"+nodeId.ToString()+" \n")
            stream.Write("CSYS, 4 \n")
            stream.Write("N \n")
            stream.Write("lastNode = ndnext(0) \n")
            stream.Write("E,lastNode,shaftNode"+geomId.ToString()+" \n")
            stream.Write("NSEL,A,,,"+nodeId.ToString()+" \n")
            stream.Write("CSYS,CSmax+1 \n")
            stream.Write("NROT,all \n")
            stream.Write("CSYS,0 \n")
            stream.Write("CP,NEXT,UX,ALL \n")
            stream.Write("CP,NEXT,UY,ALL \n")
        
        #additional displacements
        stream.Write("NSEL,,,,shaftNode"+geomId.ToString()+" \n")
        stream.Write("CSYS,CSmax+1 \n")        
        stream.Write("NROT,all \n")
        stream.Write("CSYS,0 \n")
        stream.Write("/SOLU \n")
        stream.Write("D,ALL,UZ,0 \n")
        stream.Write("D,ALL,ROTZ,0 \n")
        stream.Write("/PREP7 \n")
        
    
    stream.Write("ALLS \n")
    stream.Write("CSYS,0 \n")
    stream.Write("/SOLU \n")
    stream.Write("C*** Finished applying cartwheels \n")
    