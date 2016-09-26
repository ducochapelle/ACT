#work out a standard for adjectiveNoun or AdjectiveNoun
#TODO: named selection zodat je alleen op een subsectie kan draaien
#TODO: check if loadcase exists

from datetime import datetime
global Identifier
Identifier = {}

def ExtAPILogWriteMessage(string):
    # ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+string)
    pass

def clearLog():
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

def createLoad(analysis):
    ExtAPILogWriteMessage("createLoad...")
        
    AnalysisNames = ExtAPI.DataModel.AnalysisNames
    ExtAPILogWriteMessage(str(AnalysisNames))
    
    LoadCases = {}
    for analysisName in AnalysisNames:
        try:
            temp = ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.Result("S")
            del temp
        except:
            continue
        LoadCases[analysisName] = {}
        for set in range(1,ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.ResultSetCount+1):
            LoadCases[analysisName][set]=230
        
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
    

def doClearLog(analysis):
    clearLog()
    
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
        
        mesh = result.Analysis.MeshData
        
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
                            if result.Analysis.MeshData.MeshRegionById(edge.Id).NodeCount == 0: #Edges which are skipped with meshing
                                continue
                            if edge.Faces.Count < 2:    #If the edge doesn't have multiple faces, it's not a weld
                                continue
                            if not edge.Vertices.Count == 2:
                                continue
                            if edge.Vertices[0].Id == edge.Vertices[1].Id:
                                continue
                            NodeIdsOnEdge = mesh.MeshRegionById(edge.Id).NodeIds
                            for NodeId in NodeIdsOnEdge:     #NSLL
                                for ConnectedElementId in mesh.NodeById(NodeId).ConnectedElementIds:    #ESLN 
                                    nodesOnTheLine = 0
                                    for NodeIdeep in mesh.ElementById(ConnectedElementId).NodeIds:
                                        if NodeIdeep in NodeIdsOnEdge:
                                            nodesOnTheLine += 1
                                    if nodesOnTheLine ==1:
                                        continue
                                    if not edge.Id in [x.Id for x in EdgesByElemId[ConnectedElementId]]:
                                        EdgesByElemId[ConnectedElementId].append(edge)
                        for ElementId in mesh.MeshRegionById(face.Id).ElementIds:
                            FaceByElemId[ElementId] = face
        for elem in mesh.ElementIds:                                                    
            if not EdgesByElemId[elem].Count == 0:
                if int(mesh.ElementById(elem).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
                    WeldElements.append(elem)
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Created Geometry Dictionary...") 
        
        ExtAPILogWriteMessage(str(ExtAPI.DataModel.AnalysisNames))
        # TODO: make lc and sa read only
        lc = "LoadCases = " + str(result.Properties.GetByName("lc").Value)
        exec lc
        AnalysisNames = [analysisName for analysisName in LoadCases]
        
        EdgeByElem = result.Properties.GetByName("EdgeByElem").Value
            
        for analysisName in AnalysisNames:
            Identifier[analysisName] = {}
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            for resultSet in LoadCases[analysisName]:
                Identifier[analysisName][resultSet] = {}
                # Retrieve the stresses of the elements
                #maybe ask stress of the nodes who are actual on the line?? I think the element Ids can be easily queried...
                ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> Load Case: "+str(analysisName)+":"+str(resultSet))
                ExtAPILogWriteMessage("analysis "+str(analysisName))
                ExtAPILogWriteMessage("Get raw stresses")
                Analysis.ResultsData.CurrentResultSet = resultSet
                resultStress = Analysis.ResultsData.Result("S")
                
                allowableStress = int(result.Properties.GetByName("sa").Value)
                if allowableStress == 0:
                    allowableStress = LoadCases[analysisName][resultSet]
                
                for elemId in WeldElements:   
                    # if elemId in (92786,99611,100945): # hack
                        # continue
                    

                    ExtAPILogWriteMessage(str(datetime.now())+" >>> "+"ELEM: "+str(elemId))
                    Identifier[analysisName][resultSet][elemId]={"2xFW":{}, "2xPP":{}, "1xFW":{}, "1xPP":{}}
                    
                    plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
                    
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
                        Edge = EdgesByElemId[elemId][0]
                    else:
                        if EdgeByElem == "Envelope":
                            weldThickness = 0
                            for edge in EdgesByElemId[elemId]:
                                units = getNormalDirection(edge, FaceByElemId[elemId]) 
                                NormalStress = inProduct([sx,sy,sz],units[0])
                                # CrossStress =  inProduct([sx,sy,sz],units[2])     # HIERVOOR moet je SZY zien te vinden!!!
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
                    # Calculate the actual weld
                    units = getNormalDirection(Edge, FaceByElemId[elemId]) 
                    NormalStress = inProduct([sxMid,syMid,szMid],units[0])
                    LongStress =  inProduct([sxy,syz,sxz],units[0])  # was units[0] .... hm... does that hurt?
                    EquivalentStress = ( NormalStress**2+3*LongStress**2 )**(0.5) 
                    weldThickness = ( EquivalentStress * plateThickness ) / (allowableStress)
                    Identifier[analysisName][resultSet][elemId]["Equivalent (von-Mises) Stress"]           = abs(EquivalentStress)
                    Identifier[analysisName][resultSet][elemId]["Axial Stress"]                            = abs(NormalStress)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Stress"]                     = abs(LongStress)
                    Identifier[analysisName][resultSet][elemId]["Throat Thickness according to LasE2.mac"] = weldThickness
                    
                    NormalForce = NormalStress * plateThickness
                    LongForce = LongStress * plateThickness
                    Moment = (inProduct([sxTop,syTop,szTop],units[0]) - NormalStress) * plateThickness**2 * 6**(-1) 
                    Identifier[analysisName][resultSet][elemId]["Axial Force mm^-1"]           = abs(NormalForce)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Moment mm^-1"]   = abs(Moment)
                    Identifier[analysisName][resultSet][elemId]["Longitudinal Force mm^-1"]    = abs(LongForce)
                    
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
                    
                    #Iterate towards the minimum throat thickness with Moment uit Las Spanningen, Staal Profielen
                    
                    a = 5.0
                    for iii in range(195):
                        sigma1 = NormalForce * 1.414 / (4 * a)
                        tau1 = LongForce / (2 * a)
                        b = plateThickness + 0.67 * a * 1.414
                        sigma2 = Moment / (1.414 * a * b)
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)
                        if sigma_vm < allowableStress:
                            if abs(sigma_vm) > 1:
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            break
                        a += 1
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["2xFW"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                    
                    a = 5.0
                    for iii in range(195):
                        sigma1 = NormalForce / (2 * a)
                        tau1 = LongForce / (2 * (a / 1.414))
                        b = plateThickness - 0.67 * a
                        sigma2 = Moment / (a * b)
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)
                        if sigma_vm < allowableStress:
                            if abs(sigma_vm) > 1:
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            break
                        a += 1
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["2xPP"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                                                                                                                    
                    a = 5.0                                                                                         
                    for iii in range(195):                                                                          
                        sigma1 = NormalForce * 1.414 / (2 * a)                                                      
                        tau1 = LongForce / a                                                                        
                        sigma2 = Moment / (6**(-1) * (1.414 * a)**2)                                                
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                    
                        if sigma_vm < allowableStress:                                                              
                            if abs(sigma_vm) > 1:
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            break
                        a += 1                                                                                      
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Throat Thickness"]                                  = a
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Bending over Equivalent (von-Mises) Stress"] = bendingStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressOverVonMisses
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Equivalent (von-Mises) Stress"]              = abs(sigma_vm)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Axial Stress"]                               = abs(sigma1)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Bending Stress"]                             = abs(sigma2)
                    Identifier[analysisName][resultSet][elemId]["1xFW"]["Scaled Longitudinal Stress"]                        = abs(tau1)
                                                                                                                    
                    a = 5.0                                                                                         
                    for iii in range(195):                                                                          
                        sigma1 = NormalForce / a                                                                    
                        tau1 = LongForce / (a / 1.414)                                                              
                        sigma2 = Moment / (6**(-1) * a**2)                                                          
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                    
                        if sigma_vm < allowableStress:                                                              
                            if abs(sigma_vm) > 1:
                                bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                                longStressOverVonMisses = abs(1.732*tau1) / abs(sigma_vm)  #divide by zero
                            else:
                                bendingStressOverVonMisses = 0
                                longStressOverVonMisses = 0
                            break
                        a += 1                                                                                      
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

    Display = result.Properties.GetByName("Display").Value
    WeldType = result.Properties.GetByName("WeldType").Value
    if Display in ("Throat Thickness","Scaled Longitudinal over Equivalent (von-Mises) Stress", "Scaled Bending over Equivalent (von-Mises) Stress", "Scaled Equivalent (von-Mises) Stress", "Scaled Axial Stress", "Scaled Bending Stress", "Scaled Longitudinal Stress"):
        vals = []
        for analysisName in Identifier:
            for resultSet in Identifier[analysisName]:
                vals.append(Identifier[analysisName][resultSet][SelemId][WeldType][Display])
        return [max(vals)*1.0]
    else:
        vals = []
        for analysisName in Identifier:
            for resultSet in Identifier[analysisName]:
                vals.append(Identifier[analysisName][resultSet][SelemId][Display])
        return [max(vals)*1.0]
        # return [max([Identifier[analysisName][resultSet][SelemId][WeldType][Display] for resultSet in LoadCases[analysisName] for analysisName in AnalysisNames])]
    
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
    prop.AddOption("1xFW")
    prop.AddOption("1xPP")
    