
#work out a standard for adjectiveNoun or AdjectiveNoun
#TODO: named selection zodat je alleen op een subsectie kan draaien
#TODO: check if loadcase exists
#TODO: toch identifiers zodat er verschillende berekeningen naast elkaar kunnen draaien: S355/S690, Fatigue...
#TODO: when edge is curved; use edge nodes if they are not on top of each other.
#? When subset, allow option to enable "red" edges instead of "yellow, purple and black" ones to contain welds

from datetime import datetime
global Identifier
Identifier = {}

def ExtAPILogWriteMessage(string):
    # ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+string)
    pass

def clearLog():
    # try:
        # if CreatedReverseLookupDictionairy == 1:
            # del CreatedReverseLookupDictionairy
            # Identifier = {}
    # except:
        # pass
    
    f = open(ExtAPI.Log.LogFilename,'w')
    f.write("<html>                                        \n")
    f.write("<head>                                        \n")
    f.write('<style type="text/css">                       \n')
    f.write("  .error { color:red; }                       \n")
    f.write("</style>                                      \n")
    f.write("De tijd en datum van vandaag                  \n")
    f.write("</head>                                       \n")
    f.write("<body>                                        \n")
    f.write("<br/>                                         \n")
    f.write("<br/>                                         \n")
    f.close()

def doClearLog(analysis):
    clearLog()

def defaultLoadCases():
    LoadCases = {}
    AnalysisNames = ExtAPI.DataModel.AnalysisNames
    for analysisName in AnalysisNames:
        try:
            temp = ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.Result("S")
            del temp
        except:
            continue
        LoadCases[analysisName] = {}
        for set in range(1,ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.ResultSetCount+1):
            LoadCases[analysisName][set]=230
    return LoadCases

    
def createLoad(analysis):
    for resultObject in analysis.ResultObjects:
        if str(type(resultObject)) == "<type 'SimResult'>":
            if resultObject.Name == "MultiWeldscale":
                analysis.CreateResultObject("MultiWeldscale")
                analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("lc").Value = resultObject.Properties.GetByName("lc").Value
                # analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("sa").Value = resultObject.Properties.GetByName("sa").Value
                break
    else:
        #todo: if there are already multiweldscale objects, just add one with the same "lc" as the others
        ExtAPILogWriteMessage("createLoad...")
        LoadCases = defaultLoadCases()
        
        # todo: add names
        analysis.CreateResultObject("MultiWeldscale")
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xFW"
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("lc").Value = str(LoadCases)
        analysis.CreateResultObject("MultiWeldscale")
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xPP"
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("lc").Value = str(LoadCases)
        analysis.CreateResultObject("MultiWeldscale")
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "1xPP"
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("lc").Value = str(LoadCases)
        analysis.CreateResultObject("MultiWeldscale")
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xFW"
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("Display").Value = "Scaled Longitudinal over Equivalent (von-Mises) Stress"
        analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("lc").Value = str(LoadCases)
    
def validateLoadCases(result, property):    
    # ExtAPI.Log.WriteMessage("validate!")
    global g_result
    g_result = result
    global g_property
    g_property = property
    lc = "LoadCases = " + str(property.Value)
    try:
        exec lc
        AnalysisNames = [analysisName for analysisName in LoadCases]
        for analysisName in AnalysisNames:
            if not analysisName in ExtAPI.DataModel.AnalysisNames:
                ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" does not exist.")
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            for resultSet in LoadCases[analysisName]:
                if Analysis.ResultsData.ResultSetCount < resultSet:
                    ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" has no step "+str(resultSet)+".")
                sa = LoadCases[analysisName][resultSet]
                if sa < 10 or sa > 1000:
                    ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+":"+str(resultSet)+" has unexpected allowable: "+str(sa)+" is not between 10 and 1000. Which is already a pretty big range.")
    except:
        ExtAPI.Application.LogWarning("Invalid Load Cases.")
        property.Value = "{}"
    return
    
def getNormalDirection(edge,face):
    # returns the unitvector perpendicular to the edge and in the plane of the face
    ExtAPILogWriteMessage("start getNormalDirection ")
    ExtAPILogWriteMessage("plateNormal ")
    plateNormal = (face.Normals[0], \
                   face.Normals[1], \
                   face.Normals[2])
    ExtAPILogWriteMessage("edgeVector ")
    edgeVector = (edge.StartVertex.X - edge.EndVertex.X,
                  edge.StartVertex.Y - edge.EndVertex.Y,
                  edge.StartVertex.Z - edge.EndVertex.Z)
    ExtAPILogWriteMessage("edgeVectorLength ")
    edgeVectorLength = ( edgeVector[0]**2+edgeVector[1]**2+edgeVector[2]**2 ) ** ( 0.5 )
    ExtAPILogWriteMessage("edgeUnits ")
    edgeUnit = (edgeVector[0]/edgeVectorLength,
                edgeVector[1]/edgeVectorLength,
                edgeVector[2]/edgeVectorLength)
    ExtAPILogWriteMessage("Norm ")
    Norm = (plateNormal[1] * edgeUnit[2] - plateNormal[2] * edgeUnit[1],
            plateNormal[2] * edgeUnit[0] - plateNormal[0] * edgeUnit[2],
            plateNormal[0] * edgeUnit[1] - plateNormal[1] * edgeUnit[0])
            
    ExtAPILogWriteMessage("return getNormalDirection ")
    return (Norm,edgeUnit,plateNormal)
    
def inProduct(F,u):
    A = F[0]*u[0] + F[1]*u[1] + F[2]*u[2]
    A = A*1.0
    return A

def getValue(result,SelemId):
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+"START: GetValue of "+str(SelemId))

    ExtAPILogWriteMessage("Test existence of the fkin CreatedReverseLookupDictionairy ")
    try:
        if CreatedReverseLookupDictionairy == 1:
            ExtAPILogWriteMessage("Allready Created Reverse Lookup Dictionairy")
            pass
    except:
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Create Geometry Dictionary...")
        global CreatedReverseLookupDictionairy
        CreatedReverseLookupDictionairy = 1
        
        global FaceByElemId                                                                                 
        global EdgesByElemId 
        global WeldElements
        
        
        Mesh = result.Analysis.MeshData
        GeoData = result.Analysis.GeoData
        
        ScopeEdgeIds = []   # de shit doet niet ALLE bodies, alleen het laatste component? bij All Bodies
        for Body in [GeoData.GeoEntityById(Id) for Id in result.Properties.GetByName("Scope").Value.Ids]:
            for Face in Body.Faces:
                for Edge in Face.Edges:
                    ScopeEdgeIds.append(Edge.Id)
        ExtAPILogWriteMessage("ScopeEdgeIds: "+str(ScopeEdgeIds))            
        
        NodesOnLines = set([])                                                                              
        WeldElements = []
        FaceByElemId = {}                                                                                   
        EdgesByElemId = {}
        for elem in result.Analysis.MeshData.ElementIds:                                                    
            EdgesByElemId[elem] = []            
        for assembly in result.Analysis.GeoData.Assemblies:
            for part in assembly.Parts:
                for body in part.Bodies:
                    ExtAPILogWriteMessage("Body: "+str(body.Id))
                    if body.Suppressed == True:
                        continue
                    if not str(body.BodyType) == "GeoBodySheet":
                        continue
                    for face in body.Faces:     #ASLV
                        ExtAPILogWriteMessage("Face: "+str(face.Id))
                        for edge in face.Edges:     #LSLA
                            ExtAPILogWriteMessage("Edge: "+str(edge.Id))
                            if not edge.Id in ScopeEdgeIds:
                                ExtAPILogWriteMessage("Not in scope")
                                continue
                            if result.Analysis.MeshData.MeshRegionById(edge.Id).NodeCount == 0: #Edges which are skipped with meshing
                                ExtAPILogWriteMessage("No mesh")
                                continue
                            if edge.Faces.Count < 2:    #If the edge doesn't have multiple faces, it's not a weld
                                ExtAPILogWriteMessage("Not multiple faces")
                                continue
                            if not edge.Vertices.Count == 2:
                                ExtAPILogWriteMessage("Not two vertices")
                                continue
                            if edge.Vertices[0].Id == edge.Vertices[1].Id:
                                ExtAPILogWriteMessage("Not two distinct vertices")
                                continue
                            ExtAPILogWriteMessage("Edge seems ok. Let us proceed.")
                            NodeIdsOnEdge = Mesh.MeshRegionById(edge.Id).NodeIds
                            for NodeId in NodeIdsOnEdge:     #NSLL
                                ExtAPILogWriteMessage("NodeId: "+str(NodeId))
                                for ConnectedElementId in Mesh.NodeById(NodeId).ConnectedElementIds:    #ESLN 
                                    nodesOnTheLine = 0
                                    for NodeIdeep in Mesh.ElementById(ConnectedElementId).NodeIds:
                                        if NodeIdeep in NodeIdsOnEdge:
                                            nodesOnTheLine += 1
                                    if nodesOnTheLine ==1:
                                        continue
                                    if not edge.Id in [x.Id for x in EdgesByElemId[ConnectedElementId]]:
                                        EdgesByElemId[ConnectedElementId].append(edge)
                        for ElementId in Mesh.MeshRegionById(face.Id).ElementIds:
                            FaceByElemId[ElementId] = face
        for elem in Mesh.ElementIds:                                       
            ExtAPILogWriteMessage("Elem: "+str(elem))
            if not EdgesByElemId[elem].Count == 0:
                if int(Mesh.ElementById(elem).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
                    WeldElements.append(elem)
        ExtAPILogWriteMessage("Weld Elements: "+str(WeldElements))
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Created Geometry Dictionary...") 
        
        ExtAPILogWriteMessage(str(ExtAPI.DataModel.AnalysisNames))
        lc = "LoadCases = " + str(result.Properties.GetByName("lc").Value)
        exec lc
        AnalysisNames = [analysisName for analysisName in LoadCases]
        
        EdgeByElem = result.Properties.GetByName("EdgeByElem").Value
            
        for analysisName in AnalysisNames:
            # if not analysisName in ExtAPI.DataModel.AnalysisNames:    #checked this closer to the input
                # ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" does not exist.")
                # continue
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            Identifier[analysisName] = {}
            for resultSet in LoadCases[analysisName]:
                # if Analysis.ResultsData.ResultSetCount < resultSet:   #checked this closer to the input
                    # ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" has no step "+str(resultSet)+".")
                    # continue
                Identifier[analysisName][resultSet] = {}
                # Retrieve the stresses of the elements
                #maybe ask stress of the nodes who are actual on the line?? I think the element Ids can be easily queried...
                ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Load Case: "+str(analysisName)+":"+str(resultSet))
                ExtAPILogWriteMessage("Get stres set")
                Analysis.ResultsData.CurrentResultSet = resultSet
                resultStress = Analysis.ResultsData.Result("S")
                ExtAPILogWriteMessage("Get allowable stress")
                allowableStress = LoadCases[analysisName][resultSet]
                
                for elemId in WeldElements:   
                    ExtAPILogWriteMessage("ELEM: "+str(elemId))
                    Identifier[analysisName][resultSet][elemId]={"2xFW":{}, "2xPP":{}, "1xFW":{}, "1xPP":{}}
                    
                    plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
                    
                    ExtAPILogWriteMessage("Get raw stresses")
                    sxTotal = [abs(x) for x in resultStress.ElementValue(elemId,"X")]
                    syTotal = [abs(x) for x in resultStress.ElementValue(elemId,"Y")]
                    szTotal = [abs(x) for x in resultStress.ElementValue(elemId,"Z")]
                    sx = sum(sxTotal,0.0)/12  
                    sy = sum(syTotal,0.0)/12  
                    sz = sum(szTotal,0.0)/12  
                    sxTop = sum(sxTotal[0:4],0.0)/4  
                    syTop = sum(syTotal[0:4],0.0)/4  
                    szTop = sum(szTotal[0:4],0.0)/4  
                    sxBot = sum(sxTotal[4:8],0.0)/4
                    syBot = sum(syTotal[4:8],0.0)/4
                    szBot = sum(szTotal[4:8],0.0)/4
                    sxMid = sum(sxTotal[8:12],0.0)/4
                    syMid = sum(syTotal[8:12],0.0)/4
                    szMid = sum(szTotal[8:12],0.0)/4
                    sxy = sum([abs(x) for x in resultStress.ElementValue(elemId,"XY")],0.0)/12 
                    syz = sum([abs(x) for x in resultStress.ElementValue(elemId,"YZ")],0.0)/12 
                    sxz = sum([abs(x) for x in resultStress.ElementValue(elemId,"XZ")],0.0)/12 
                    
                    # Determine the edge to calculate with; for elements can be with multiple edges
                    # Contains duplicate calculations to be sure the worst-case edge is taken into account
                    if EdgesByElemId[elemId].Count == 1:
                        ExtAPILogWriteMessage("Only one edge")
                        Edge = EdgesByElemId[elemId][0]
                    else:
                        ExtAPILogWriteMessage("Determine governing edge")
                        if EdgeByElem == "Envelope":
                            weldThickness = 0
                            for edge in EdgesByElemId[elemId]:
                                units = getNormalDirection(edge, FaceByElemId[elemId]) 
                                NormalStress = inProduct([sx,sy,sz],units[0])
                                LongStress =  inProduct([sxy,syz,sxz],units[1]) 
                                newWeldThickness = ( ( NormalStress**2+3*LongStress**2 )**(0.5) * plateThickness ) / (allowableStress)
                                if newWeldThickness >= weldThickness:
                                    EdgeWithMaxWeldThickness = edge
                                    weldThickness = newWeldThickness
                            Edge = EdgeWithMaxWeldThickness
                        elif EdgeByElem == "High Edge Ids":
                            Edge = max([x.Id for x in EdgesByElemId[elemId]])
                        elif EdgeByElem == "Low Edge Ids":
                            Edge = min([x.Id for x in EdgesByElemId[elemId]])
                    ExtAPILogWriteMessage("Calculate the actual weld")
                    ExtAPILogWriteMessage("1")
                    # Calculate the actual weld
                    units = getNormalDirection(Edge, FaceByElemId[elemId])  # don't use edge but edge nodes? Check if edge is GeoCurveLine or GeoCurveCircle
                    ExtAPILogWriteMessage("2")
                    NormalStress = inProduct([sxMid,syMid,szMid],units[0])
                    ExtAPILogWriteMessage("3")
                    LongStress =  inProduct([sxy,syz,sxz],units[0])  
                    ExtAPILogWriteMessage("4")
                    EquivalentStress = ( NormalStress**2+3*LongStress**2 )**(0.5) 
                    ExtAPILogWriteMessage("5")
                    weldThickness = ( EquivalentStress * plateThickness ) / (allowableStress)
                    ExtAPILogWriteMessage("6")
                    Identifier[analysisName][resultSet][elemId]["Equivalent (von-Mises) Stress"]           = abs(EquivalentStress)
                    Identifier[analysisName][resultSet][elemId]["Axial Stress"]                            = abs(NormalStress)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Stress"]                     = abs(LongStress)
                    Identifier[analysisName][resultSet][elemId]["Throat Thickness according to LasE2.mac"] = weldThickness
                    ExtAPILogWriteMessage("7")
                    
                    NormalForce = abs( NormalStress * plateThickness )
                    ExtAPILogWriteMessage("8")
                    LongForce = abs( LongStress * plateThickness )
                    ExtAPILogWriteMessage("9")
                    Moment = abs( (inProduct([sxTop,syTop,szTop],units[0]) - NormalStress) * plateThickness**2 * 6**(-1) )
                    ExtAPILogWriteMessage("A")
                    Identifier[analysisName][resultSet][elemId]["Axial Force mm^-1"]           = abs(NormalForce)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Moment mm^-1"]   = abs(Moment)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Force mm^-1"]    = abs(LongForce)
                    ExtAPILogWriteMessage("B")
                    
                    # Get unit vectors for debuggin purpose
                    ExtAPILogWriteMessage("Set normal units")
                    unitNormalScalar = round(abs(9*units[0][0]))*100 \
                                     + round(abs(9*units[0][1]))*10 \
                                     + round(abs(9*units[0][2]))
                    unitCrossScalar  = round(abs(9*units[2][0]))*100 \
                                     + round(abs(9*units[2][1]))*10 \
                                     + round(abs(9*units[2][2]))
                    unitLongScalar   = round(abs(9*units[1][0]))*100 \
                                     + round(abs(9*units[1][1]))*10 \
                                     + round(abs(9*units[1][2]))
                    Identifier[analysisName][resultSet][elemId]["Axial Unit Vector"]           = unitNormalScalar
                    Identifier[analysisName][resultSet][elemId]["Cross Unit Vector"]           = unitCrossScalar
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Unit Vector"]    = unitLongScalar
                    ExtAPILogWriteMessage("Setted normal units")
                    
                    #Iterate towards the minimum throat thickness with Moment uit Las Spanningen, Staal Profielen
                    ExtAPILogWriteMessage("2xFW")
                    a = 5.0
                    for iii in range(195):
                        ExtAPILogWriteMessage("a = "+str(a))
                        sigma1 = NormalForce * 1.414 / (4 * a)
                        tau1 = LongForce / (2 * a)
                        b = plateThickness + 0.67 * a * 1.414
                        sigma2 = Moment / (1.414 * a * b)
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 + 3 * (sigma1 + sigma2)**2 )**(0.5)
                        ExtAPILogWriteMessage("1.1")
                        if sigma_vm < allowableStress:
                            if abs(sigma_vm) > 1:
                                ExtAPILogWriteMessage("1.2")
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                ExtAPILogWriteMessage("1.3")
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                ExtAPILogWriteMessage("1.4")
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            ExtAPILogWriteMessage("1.5")
                            break
                        ExtAPILogWriteMessage("1.6")
                        a += 1
                    ExtAPILogWriteMessage("1.7")
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Axial Stress"]                               = abs(sigma1 * 2.732)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Bending Stress"]                             = abs(sigma2 * 2.732)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                    
                    ExtAPILogWriteMessage("2xPP")
                    a = 5.0
                    for iii in range(195):
                        ExtAPILogWriteMessage("a = "+str(a))
                        sigma1 = NormalForce / (2 * a)
                        tau1 = LongForce / (2 * (a / 1.414))
                        b = plateThickness - a
                        if b == 0:
                            continue
                        sigma2 = Moment / (a * b)   # divide by zero
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)
                        ExtAPILogWriteMessage("2.1")
                        if sigma_vm < allowableStress:
                            if abs(sigma_vm) > 1:
                                ExtAPILogWriteMessage("2.2")
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                ExtAPILogWriteMessage("2.3")
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                ExtAPILogWriteMessage("2.4")
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            ExtAPILogWriteMessage("2.5")
                            break
                        ExtAPILogWriteMessage("2.6")
                        a += 1
                    ExtAPILogWriteMessage("2.7")
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                                                                                                                    
                    ExtAPILogWriteMessage("1xFW")
                    a = 5.0                                                                                         
                    for iii in range(195):                                                                          
                        ExtAPILogWriteMessage("a = "+str(a))
                        sigma1 = NormalForce * 1.414 / (2 * a)                                                      
                        tau1 = LongForce / a                                                                        
                        sigma2 = Moment / (6**(-1) * (1.414 * a)**2)                                                
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)
                        ExtAPILogWriteMessage("3.1")
                        if sigma_vm < allowableStress:                                                              
                            ExtAPILogWriteMessage("3.2")
                            if abs(sigma_vm) > 1:
                                ExtAPILogWriteMessage("3.3")
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                ExtAPILogWriteMessage("3.4")
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            ExtAPILogWriteMessage("3.5")
                            break
                        ExtAPILogWriteMessage("3.6")
                        a += 1                                                                                      
                    ExtAPILogWriteMessage("3.7")
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                                                                                                                    
                    ExtAPILogWriteMessage("1xPP")
                    a = 5.0                                                                                         
                    for iii in range(195):                                                                          
                        ExtAPILogWriteMessage("a = "+str(a))
                        sigma1 = NormalForce / a                                                                    
                        tau1 = LongForce / (a / 1.414)                                                              
                        sigma2 = Moment / (6**(-1) * a**2)                                                          
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)
                        ExtAPILogWriteMessage("4.1")
                        if sigma_vm < allowableStress:                                                              
                            ExtAPILogWriteMessage("4.2")
                            if abs(sigma_vm) > 1:
                                ExtAPILogWriteMessage("4.3")
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                ExtAPILogWriteMessage("4.4")
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            ExtAPILogWriteMessage("4.5")
                            break
                        ExtAPILogWriteMessage("4.6")
                        a += 1                                                                                      
                    ExtAPILogWriteMessage("4.7")
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["1xPP"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                    ExtAPILogWriteMessage(str(Identifier[analysisName][resultSet][elemId]))
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Created Results Dictionary...")
    if SelemId not in WeldElements:
        ExtAPILogWriteMessage("CONTINUE: elemId not in WeldElements")
        return []
    ExtAPILogWriteMessage("Get Display")
    Display = result.Properties.GetByName("Display").Value
    ExtAPILogWriteMessage("Get WeldType")
    WeldType = result.Properties.GetByName("WeldType").Value
    if Display in ("Throat Thickness","Scaled Longitudinal over Equivalent (von-Mises) Stress", "Scaled Bending over Equivalent (von-Mises) Stress", "Scaled Equivalent (von-Mises) Stress", "Scaled Axial Stress", "Scaled Bending Stress", "Scaled Longitudinal Stress"):
        ExtAPILogWriteMessage("Create first max list")
        vals = []
        for analysisName in Identifier:
            for resultSet in Identifier[analysisName]:
                vals.append(Identifier[analysisName][resultSet][SelemId][WeldType][Display])
        ExtAPILogWriteMessage("Created first max list")
        return [max(vals)*1.0]
    else:
        ExtAPILogWriteMessage("Create second max list")
        vals = []
        for analysisName in Identifier:
            for resultSet in Identifier[analysisName]:
                vals.append(Identifier[analysisName][resultSet][SelemId][Display])
        ExtAPILogWriteMessage("Created second max list")
        return [max(vals)*1.0]
    
def selectDisplay(load,prop):
    prop.ClearOptions()
    prop.AddOption("Throat Thickness")
    prop.AddOption(" ")
    prop.AddOption("Scaled Equivalent (von-Mises) Stress")
    prop.AddOption("Scaled Bending over Equivalent (von-Mises) Stress")
    prop.AddOption("Scaled Longitudinal over Equivalent (von-Mises) Stress")
    prop.AddOption("Scaled Axial Stress")
    prop.AddOption("Scaled Bending Stress")
    prop.AddOption("Scaled Longitudinal Stress")
    prop.AddOption("Axial Force mm^-1")
    prop.AddOption("Longitudinal Moment mm^-1")
    prop.AddOption("Longitudinal Force mm^-1")
    prop.AddOption("Equivalent (von-Mises) Stress")
    prop.AddOption("Axial Stress")
    prop.AddOption("Longitudinal Stress")
    prop.AddOption("Axial Unit Vector")
    prop.AddOption("Cross Unit Vector")
    prop.AddOption("Longitudinal Unit Vector")
    prop.AddOption("Throat Thickness according to LasE2.mac")
    

def selectEdgeByElem(load,prop):
    prop.ClearOptions()
    prop.AddOption("Envelope")
    prop.AddOption("High Edge Ids")
    prop.AddOption("Low Edge Ids")

def selectWeldType(load,prop):  
    prop.ClearOptions()
    prop.AddOption("2xFW")
    prop.AddOption("2xPP")
    # prop.AddOption("1xFW")
    prop.AddOption("1xPP")
    
def validateScope(result, prop):
    if not result.Analysis.GeoData.GeoEntityById(prop.Value.Ids[0]).Type == GeoCellTypeEnum.GeoBody:
        prop.Value.Ids = []
    
    