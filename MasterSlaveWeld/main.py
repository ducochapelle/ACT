#work out a standard for adjectiveNoun or AdjectiveNoun

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
    analysis.CreateResultObject("MultiWeldscale")
    analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xFW"
    analysis.CreateResultObject("MultiWeldscale")
    analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xPP"
    analysis.CreateResultObject("MultiWeldscale")
    analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "1xPP"
    analysis.CreateResultObject("MultiWeldscale")
    analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("WeldType").Value = "2xFW"
    analysis.ResultObjects[analysis.ResultObjects.Count-1].Properties.GetByName("Display").Value = "Scaled Longitudinal over Equivalent (von-Mises) Stress"

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
            
    # edgeStartVertex = edge.StartVertex
    # edgeEndVertex = edge.EndVertex
    
    # plateNormal = face.Normals[0:3]
    
    # edgeVector = (edgeStartVertex.X - edgeEndVertex.X,
                  # edgeStartVertex.Y - edgeEndVertex.Y,
                  # edgeStartVertex.Z - edgeEndVertex.Z)
    # edgeVectorLength = ( edgeVector[0]**2+edgeVector[1]**2+edgeVector[2]**2 ) ** ( 0.5 )
    # edgeUnit = (edgeVector[0]/edgeVectorLength,
                # edgeVector[1]/edgeVectorLength,
                # edgeVector[2]/edgeVectorLength)
    # Norm = (plateNormal[1] * edgeUnit[2] - plateNormal[2] * edgeUnit[1],
            # plateNormal[2] * edgeUnit[0] - plateNormal[0] * edgeUnit[2],
            # plateNormal[0] * edgeUnit[1] - plateNormal[1] * edgeUnit[0])
            
    # return {"Normal":Norm, "Long":edgeUnit, "Cross":plateNormal}      # tested for performance, dict vs tuple, didnt matter. tuple vs array doesnt matter either.
    ExtAPILogWriteMessage("return getNormalDirection ")
    return (Norm,edgeUnit,plateNormal)
    
def inProduct(F,u):
    A = F[0]*u[0] + F[1]*u[1] + F[2]*u[2]
    A = A*1.0
    return A

def getValue(result,SelemId):
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+"START: GetValue of "+str(SelemId))
    


    ExtAPILogWriteMessage("Test existence of CreatedReverseLookupDictionairy ")
    try:
        if CreatedReverseLookupDictionairy == 1:
            ExtAPILogWriteMessage("Allready Created Reverse Lookup Dictionairy")
            pass
    except:
        # push this whole except in a seperate button on the toolbar?
        ExtAPILogWriteMessage("Create Reverse Lookup Dictionairy ")
        global CreatedReverseLookupDictionairy
        CreatedReverseLookupDictionairy = 1
        
        global FaceByElemId                                                                                 
        global EdgesByElemId 
        global WeldElements
        
        NodesOnLines = set([])                                                                              
        WeldElements = []
        FaceByElemId = {}                                                                                   
        EdgesByElemId = {}
        for elem in result.Analysis.MeshData.ElementIds:                                                    
            EdgesByElemId[elem] = []            
        for assembly in result.Analysis.GeoData.Assemblies:
            for part in assembly.Parts:
                for body in part.Bodies:
                    if body.Suppressed == True:
                        continue
                    for face in body.Faces:     #ASLV
                        for edge in face.Edges:     #LSLA
                            if result.Analysis.MeshData.MeshRegionById(edge.Id).NodeCount == 0: #Edges which are skipped with meshing
                                continue
                            elif edge.Faces.Count < 2:    #If the edge doesn't have multiple faces, it's not a weld
                                continue
                            elif edge.Vertices.Count == 0:
                                continue
                            for NodeId in result.Analysis.MeshData.MeshRegionById(edge.Id).NodeIds:     #NSLL
                                for ConnectedElementId in result.Analysis.MeshData.NodeById(NodeId).ConnectedElementIds:    #ESLN 
                                    if not edge.Id in [x.Id for x in EdgesByElemId[ConnectedElementId]]:
                                        EdgesByElemId[ConnectedElementId].append(edge)
                        for ElementId in result.Analysis.MeshData.MeshRegionById(face.Id).ElementIds:
                            FaceByElemId[ElementId] = face
        ExtAPILogWriteMessage("Created Reverse Lookup Dictionairy ") 
        for elem in result.Analysis.MeshData.ElementIds:                                                    
            if not EdgesByElemId[elem].Count == 0:
                if int(result.Analysis.MeshData.ElementById(SelemId).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
                    WeldElements.append(elem)

        # ExtAPILogWriteMessage(str(ExtAPI.DataModel.AnalysisNames))

        global AnalysisNames
        global LoadCases

        
        AnalysisNames = ExtAPI.DataModel.AnalysisNames
        ExtAPILogWriteMessage(str(AnalysisNames))
        
        LoadCases = {}
        for analysisName in AnalysisNames:
            LoadCases[analysisName] = range(1,ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.ResultSetCount+1)
            
        for analysisName in AnalysisNames:
            Identifier[analysisName] = {}
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            for resultSet in LoadCases[analysisName]:
                Identifier[analysisName][resultSet] = {}
                # Retrieve the stresses of the elements
                #maybe ask stress of the nodes who are actual on the line?? I think the element Ids can be easily queried...
                
                ExtAPILogWriteMessage("analysis "+str(analysisName))
                ExtAPILogWriteMessage("Get raw stresses")
                Analysis.ResultsData.CurrentResultSet = resultSet
                resultStress = Analysis.ResultsData.Result("S")
                for elemId in ExtAPI.DataModel.MeshDataByName("Global").ElementIds:
                    
                    ExtAPILogWriteMessage("ELEM: "+str(elemId)) 
                    Identifier[analysisName][resultSet][elemId]={"2xFW":{}, "2xPP":{}, "1xFW":{}, "1xPP":{}}

                    if elemId not in WeldElements:
                        ExtAPILogWriteMessage("CONTINUE: elemId not in WeldElements")
                        continue
                    # if EdgesByElemId[elemId].Count == 0:
                        # ExtAPILogWriteMessage("CONTINUE: "+str(elemId)+" Element is not a weld") 
                        # continue
                    # if not int(result.Analysis.MeshData.ElementById(SelemId).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
                        # ExtAPILogWriteMessage("CONTINUE: Element is not a shell") 
                        # continue

                    plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
                    allowableStress = int(result.Properties.GetByName("sa").Value)

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
                    if result.Properties.GetByName("EdgeByElem").Value == "Envelope":
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
                    elif result.Properties.GetByName("EdgeByElem").Value == "High Edge Ids":
                        Edge = max([x.Id for x in EdgesByElemId[elemId]])
                    elif result.Properties.GetByName("EdgeByElem").Value == "Low Edge Ids":
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

    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> 01")
            
    if SelemId not in WeldElements:
        ExtAPILogWriteMessage("CONTINUE: elemId not in WeldElements")
        return []
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> 02")
    # if not int(result.Analysis.MeshData.ElementById(SelemId).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
        # ExtAPILogWriteMessage("RETURN: Element is not a shell") 
        # return []
    # if EdgesByElemId[SelemId].Count == 0:
        # ExtAPILogWriteMessage("RETURN: Element is not a weld") 
        # return []
    Display = result.Properties.GetByName("Display").Value
    WeldType = result.Properties.GetByName("WeldType").Value
    if Display in ("Throat Thickness","Scaled Longitudinal over Equivalent (von-Mises) Stress", "Scaled Bending over Equivalent (von-Mises) Stress", "Scaled Equivalent (von-Mises) Stress", "Scaled Axial Stress", "Scaled Bending Stress", "Scaled Longitudinal Stress"):
        vals = []
        for analysisName in AnalysisNames:
            ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> 04 "+str(analysisName))
            for resultSet in LoadCases[analysisName]:
                ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> 05 "+str(resultSet))
                vals.append(Identifier[analysisName][resultSet][SelemId][WeldType][Display])
                ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> 06"+str(vals))
        return [max(vals)]
    else:
        vals = []
        for analysisName in AnalysisNames:
            for resultSet in LoadCases[analysisName]:
                vals.append(Identifier[analysisName][resultSet][SelemId][Display])
        return [max(vals)]
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
    