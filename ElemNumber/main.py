#add node selection?
if ExtAPI.Context == "Mechanical":
    ctxCartwheel = None
    LocationByNodeId = None

def createElemNumber(analysis):
    analysis.CreateResultObject("ElemNumber")

def getElemNumber(result,elemId):
    e = int(result.Properties.GetByName("ElemNumber").Value)
    if e == 0:
        return [elemId*1.0]
    # elif e == elemId:
        # return [elemId*1.0]
    return []

def showElement(result):
        
    elemId = int(result.Properties.GetByName("ElemNumber").Value)
    if elemId == 0:
        return
    mesh = ExtAPI.DataModel.MeshDataByName("Global")
    if not elemId in mesh.ElementIds:
        ExtAPI.Application.LogWarning("Element "+str(elemId)+" does not exist.")
        return
    NodeIds = mesh.ElementById(elemId).CornerNodeIds
    global ctxCartwheel
    if ctxCartwheel!=None:
        ctxCartwheel.Visible = False
    
    LocationByNodeId = {}
    for NodeId in NodeIds:
        LocationByNodeId[NodeId] = ( mesh.NodeById(NodeId).X,mesh.NodeById(NodeId).Y,mesh.NodeById(NodeId).Z )
    
    #draw nice counters depending on element type
    
    pts = System.Array.CreateInstance(float,6)   
    ctxCartwheel = ExtAPI.Graphics.CreateAndOpenDraw3DContext()
    ctxCartwheel.Color = 0xFF0000
    ctxCartwheel.LineWeight = 2
    
    if str(result.Properties.GetByName("Style").Value) == "Fancy":
        # kBeam3 (3), kBeam4 (4), kHex20 (12), kHex8 (11), kLine2 (1), kLine3 (2), kPoint0 (0), kPyramid13 (16), kPyramid5 (15), kQuad4 (7), kQuad8 (8), kTet10 (10), kTet4 (9), kTri3 (5), kTri6 (6), kUnknown (-1), kWedge15 (14), kWedge6 (13)

        elemType = int(ExtAPI.DataModel.MeshDataByName("Global").ElementById(elemId).Type)
        Connector = None
        if elemType in (12, 20):
            Connector = ((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))
        elif elemType in (9, 10):
            Connector = ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
        elif elemType in (7, 8):
            Connector = ((0,1),(1,2),(2,3),(3,0))
        elif elemType in (5, 6):
            Connector = ((0,1),(1,2),(2,0))
        if Connector:
            for tuple in Connector:
                pts[0]=LocationByNodeId[NodeIds[tuple[0]]][0]
                pts[1]=LocationByNodeId[NodeIds[tuple[0]]][1]
                pts[2]=LocationByNodeId[NodeIds[tuple[0]]][2]
                pts[3]=LocationByNodeId[NodeIds[tuple[1]]][0]
                pts[4]=LocationByNodeId[NodeIds[tuple[1]]][1]
                pts[5]=LocationByNodeId[NodeIds[tuple[1]]][2]
                ctxCartwheel.DrawPolyline(pts)
        
    else:
        for NodeId in NodeIds:
            pts[3] = pts[0]
            pts[4] = pts[1]
            pts[5] = pts[2]
            pts[0]=LocationByNodeId[NodeId][0]
            pts[1]=LocationByNodeId[NodeId][1]
            pts[2]=LocationByNodeId[NodeId][2]
            ctxCartwheel.DrawPolyline(pts)

    ctxCartwheel.Close()
    ctxCartwheel.Visible = True  
    return        

def hideElement(a):
    global ctxCartwheel
    if ctxCartwheel!=None:
        ctxCartwheel.Visible = False
    ctxCartwheel = None
    pass

def validateElement(a,b):
    hideElement(a)
    showElement(a)
    return
    
def selectStyle(analysis,prop):
    prop.ClearOptions()
    prop.AddOption("Fancy")
    prop.AddOption("Shabby")
