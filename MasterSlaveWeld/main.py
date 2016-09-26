from datetime import datetime
global Identifier
Identifier = {}

def ExtAPILogWriteMessage(string):
    ExtAPI.Log.WriteMessage(string)
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

def doClearLog(analysis):
    clearLog()
    
def getNormalDirection(edge,face):
    # returns the unitvector perpendicular to the edge and in the plane of the face
    ExtAPILogWriteMessage(str(datetime.now())+" >>> start getNormalDirection ")
    ExtAPILogWriteMessage(str(datetime.now())+" >>> plateNormal ")
    plateNormal = (face.Normals[0], \
                   face.Normals[1], \
                   face.Normals[2])
    ExtAPILogWriteMessage(str(datetime.now())+" >>> edgeVector ")
    edgeVector = (edge.StartVertex.X - edge.EndVertex.X,
                  edge.StartVertex.Y - edge.EndVertex.Y,
                  edge.StartVertex.Z - edge.EndVertex.Z)
    ExtAPILogWriteMessage(str(datetime.now())+" >>> edgeVectorLength ")
    edgeVectorLength = ( edgeVector[0]**2+edgeVector[1]**2+edgeVector[2]**2 ) ** ( 0.5 )
    ExtAPILogWriteMessage(str(datetime.now())+" >>> edgeUnits ")
    edgeUnit = (edgeVector[0]/edgeVectorLength,
                edgeVector[1]/edgeVectorLength,
                edgeVector[2]/edgeVectorLength)
    ExtAPILogWriteMessage(str(datetime.now())+" >>> Norm ")
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
    ExtAPILogWriteMessage(str(datetime.now())+" >>> return getNormalDirection ")
    return (Norm,edgeUnit,plateNormal)
    
def inProduct(F,u):
    A = F[0]*u[0] + F[1]*u[1] + F[2]*u[2]
    A = A*1.0
    return A

def getValue(result,SelemId):
    ExtAPILogWriteMessage(str(datetime.now())+" >>> START: GetValue of "+str(SelemId))
    
    ExtAPILogWriteMessage(str(datetime.now())+" >>> Test existence of CreatedReverseLookupDictionairy ")
    try:
        if CreatedReverseLookupDictionairy == 1:
            ExtAPILogWriteMessage(str(datetime.now())+" >>> Allready Created Reverse Lookup Dictionairy")
            pass
    except:
        ExtAPILogWriteMessage(str(datetime.now())+" >>> Create Reverse Lookup Dictionairy ")
        global CreatedReverseLookupDictionairy
        CreatedReverseLookupDictionairy = 1
        
        
        global FaceByElemId                                                                                 
        global EdgesByElemId 
        
        NodesOnLines = set([])                                                                              
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
        ExtAPILogWriteMessage(str(datetime.now())+" >>> Created Reverse Lookup Dictionairy ") 
        ExtAPILogWriteMessage(str(datetime.now())+str(ExtAPI.DataModel.AnalysisNames))

        for analysisName in ExtAPI.DataModel.AnalysisNames:
            Identifier[analysisName] = {}
        
        for analysisName in ExtAPI.DataModel.AnalysisNames:
            # Retrieve the stresses of the elements
            #maybe ask stress of the nodes who are actual on the line?? I think the element Ids can be easily queried...
            
            ExtAPILogWriteMessage(str(datetime.now())+" >>> analysis "+str(analysisName))
            ExtAPILogWriteMessage(str(datetime.now())+" >>> Get raw stresses")
            resultStress = ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.Result("S")
            for elemId in result.Analysis.MeshData.ElementIds:
                
                ExtAPILogWriteMessage(str(datetime.now())+" >>> ELEM: "+str(elemId)) 
                vals = []
                
                if EdgesByElemId[elemId].Count == 0:
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> RETURN: "+str(elemId)+" Element is not a weld") 
                    continue

                ExtAPILogWriteMessage(str(datetime.now())+" >>> Get platethickness")
                plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
                ExtAPILogWriteMessage(str(datetime.now())+" >>> "+str(plateThickness))
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Get Allowablestress")
                allowableStress = int(result.Properties.GetByName("sa").Value)

                sxTotal = [abs(x) for x in resultStress.ElementValue(elemId,"X")]
                syTotal = [abs(x) for x in resultStress.ElementValue(elemId,"Y")]
                szTotal = [abs(x) for x in resultStress.ElementValue(elemId,"Z")]
                ExtAPILogWriteMessage(str(datetime.now())+" "+str((sxTotal, syTotal, szTotal)))
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
                ExtAPILogWriteMessage(str(datetime.now())+" "+str((sxy, syz, sxz)))

                
                # Determine the edge to calculate with; for elements can be with multiple edges
                # Contains duplicate calculations to be sure the worst-case edge is taken into account
                ExtAPILogWriteMessage(str(datetime.now())+" >>> start the whole extra edges misery")
                if result.Properties.GetByName("EdgeByElem").Value == "Envelope":
                    weldThickness = 0
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> Edges: "+str([x.Id for x in EdgesByElemId[elemId]]))
                    for edge in EdgesByElemId[elemId]:
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> Run the edgeId tester with edge "+str(edge.Id))
                        units = getNormalDirection(edge, FaceByElemId[elemId]) 
                        NormalStress = inProduct([sx,sy,sz],units[0])
                        # CrossStress =  inProduct([sx,sy,sz],units[2])     # HIERVOOR moet je SZY zien te vinden!!!
                        LongStress =  inProduct([sxy,syz,sxz],units[1]) 
                        newWeldThickness = ( ( NormalStress**2+3*LongStress**2 )**(0.5) * plateThickness ) / (allowableStress)
                        if newWeldThickness >= weldThickness:
                            EdgeIdWithMaxWeldThickness = edge.Id
                            weldThickness = newWeldThickness
                    EdgeId = EdgeIdWithMaxWeldThickness
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> Ran the edgeId tester witht the winner edge: "+str(EdgeId))
                elif result.Properties.GetByName("EdgeByElem").Value == "High Edge Ids":
                    EdgeId = max([x.Id for x in EdgesByElemId[elemId]])
                elif result.Properties.GetByName("EdgeByElem").Value == "Low Edge Ids":
                    EdgeId = min([x.Id for x in EdgesByElemId[elemId]])
                ExtAPILogWriteMessage(str(datetime.now())+" >>> end the whole extra edges misery")
                # Calculate the actual weld
                ExtAPILogWriteMessage(str(datetime.now())+" >>> create dummy var")
                aap = edge
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate units")
                units = getNormalDirection(aap, FaceByElemId[elemId]) 
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate normal stress")
                NormalStress = inProduct([sxMid,syMid,szMid],units[0])
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate long stress")
                LongStress =  inProduct([sxy,syz,sxz],units[0])  # was units[0] .... hm... does that hurt?
                ExtAPILogWriteMessage(str(datetime.now())+" >>> "+str((NormalStress,LongStress)))
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate equi stress")
                EquivalentStress = ( NormalStress**2+3*LongStress**2 )**(0.5) 
                ExtAPILogWriteMessage(str(datetime.now())+" >>> "+str(EquivalentStress))
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate weld thickness raw stress")
                weldThickness = ( EquivalentStress * plateThickness ) / (allowableStress)
             
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate normal and long force")
                NormalForce = NormalStress * plateThickness
                LongForce = LongStress * plateThickness
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate moment")
                Moment = (inProduct([sxTop,syTop,szTop],units[0]) - NormalStress) * plateThickness**2 * 6**(-1) 
                ExtAPILogWriteMessage(str(datetime.now())+" >>> "+str((NormalForce, LongForce, Moment)))

                # Get unit vectors for debuggin purpose
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Set normal units")
                unitNormalScalar = round(abs(9*units[0][0]))*100 \
                                 + round(abs(9*units[0][1]))*10 \
                                 + round(abs(9*units[0][2]))
                unitCrossScalar  = round(abs(9*units[2][0]))*100 \
                                 + round(abs(9*units[2][1]))*10 \
                                 + round(abs(9*units[2][2]))
                unitLongScalar   = round(abs(9*units[1][0]))*100 \
                                 + round(abs(9*units[1][1]))*10 \
                                 + round(abs(9*units[1][2]))
                
                
                #Iterate towards the minimum throat thickness with Moment
                if result.Properties.GetByName("WeldType").Value == "2xFW":
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculate 2xFW")
                    a = 5.0
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> 01 02 03")
                    for iii in range(200):
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 01")
                        sigma1 = NormalForce * 1.414 / (4 * a)                                                    # uit Las Spanningen, Staal Profielen
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 02")
                        tau1 = LongForce / (2 * a)                                                                     # uit Las Spanningen, Staal Profielen
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 03")
                        b = plateThickness + 0.67 * a * 1.414                                                   # uit Las Spanningen, Staal Profielen
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 04")
                        sigma2 = Moment / (1.414 * a * b)                                                      # uit Las Spanningen, Staal Profielen
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 05")
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                     # uit Las Spanningen, Staal Profielen
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 06")
                        if sigma_vm < allowableStress:
                            ExtAPILogWriteMessage(str(datetime.now())+" >>> 07")
                            # bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                            ExtAPILogWriteMessage(str(datetime.now())+" >>> 08")
                            weldThicknessWithMoment = a
                            ExtAPILogWriteMessage(str(datetime.now())+" >>> 09")
                            break
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 10")
                        a += 1
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 11")
                    else:
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 12")
                        weldThicknessWithMoment = 200.0
                        ExtAPILogWriteMessage(str(datetime.now())+" >>> 13")
                    ExtAPILogWriteMessage(str(datetime.now())+" >>> Calculated 2xFW")
                elif result.Properties.GetByName("WeldType").Value == "2xPP":
                    a = 5.0
                    for iii in range(200):
                        sigma1 = NormalForce / (2 * a)                                                    
                        tau1 = LongForce / (2 * (a / 1.414))                                                                    
                        b = plateThickness - 0.67 * a                                               
                        sigma2 = Moment / (a * b)                                                   
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                
                        if sigma_vm < allowableStress:
                            # bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                            weldThicknessWithMoment = a
                            break
                        a += 1
                    else:
                        weldThicknessWithMoment = 200.0
                elif result.Properties.GetByName("WeldType").Value == "1xFW":
                    a = 5.0
                    for iii in range(200):
                        sigma1 = NormalForce * 1.414 / (2 * a)                                                    # uit Las Spanningen, Staal Profielen
                        tau1 = LongForce / a                                                                     # uit Las Spanningen, Staal Profielen
                        sigma2 = Moment / (6**(-1) * (1.414 * a)**2)                                                      # uit Las Spanningen, Staal Profielen
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                     # uit Las Spanningen, Staal Profielen
                        if sigma_vm < allowableStress:
                            # bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                            weldThicknessWithMoment = a
                            break
                        a += 1
                    else:
                        weldThicknessWithMoment = 200.0
                elif result.Properties.GetByName("WeldType").Value == "1xPP":
                    a = 5.0
                    for iii in range(200):
                        sigma1 = NormalForce / a                                                    
                        tau1 = LongForce / (a / 1.414)                                                                    
                        sigma2 = Moment / (6**(-1) * a**2)                                                   
                        sigma_vm = ( (sigma1 + sigma2)**2 + 3 * tau1**2 )**(0.5)                                
                        if sigma_vm < allowableStress:
                            # bendingStressOverVonMisses = abs(sigma2) / abs(sigma_vm)  #divide by zero
                            weldThicknessWithMoment = a
                            break
                        a += 1
                    else:
                        weldThicknessWithMoment = 200.0
                            
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Create vals")
                # Output                   
                if result.Properties.GetByName("Display").Value == "Throat Thickness":
                    vals.Add(abs(weldThicknessWithMoment))
                elif result.Properties.GetByName("Display").Value == "Throat Thickness according to LasE2.mac":
                    vals.Add(weldThickness)   
                elif result.Properties.GetByName("Display").Value == "Scaled Bending over Equivalent (von-Mises) Stress":
                    vals.Add(bendingStressOverVonMisses)   
                elif result.Properties.GetByName("Display").Value == "Scaled Equivalent (von-Mises) Stress":
                    vals.Add(abs(sigma_vm))
                elif result.Properties.GetByName("Display").Value == "Scaled Axial Stress":
                    vals.Add(abs(sigma1))
                elif result.Properties.GetByName("Display").Value == "Scaled Bending Stress":
                    vals.Add(abs(sigma2))
                elif result.Properties.GetByName("Display").Value == "Scaled Longitudinal Stress":
                    vals.Add(abs(tau1))
                elif result.Properties.GetByName("Display").Value == "Axial Force mm^-1":
                    vals.Add(abs(NormalForce))
                elif result.Properties.GetByName("Display").Value == "Longitudinal Moment mm^-1":
                    vals.Add(abs(Moment))
                elif result.Properties.GetByName("Display").Value == "Longitudinal Force mm^-1":
                    vals.Add(abs(LongForce))
                elif result.Properties.GetByName("Display").Value == "Equivalent (von-Mises) Stress":
                    vals.Add(abs(EquivalentStress))
                elif result.Properties.GetByName("Display").Value == "Axial Stress":
                    vals.Add(abs(NormalStress))
                elif result.Properties.GetByName("Display").Value == "Longitudinal Stress":
                    vals.Add(abs(LongStress))
                elif result.Properties.GetByName("Display").Value == "Axial Unit Vector":
                    vals.Add(unitNormalScalar)
                elif result.Properties.GetByName("Display").Value == "Cross Unit Vector":
                    vals.Add(unitCrossScalar)
                elif result.Properties.GetByName("Display").Value == "Longitudinal Unit Vector":
                    vals.Add(unitLongScalar)
                    
                # ---- Identifier Slave
                # Make identifier name a combination of Environment name and weld type by default
                # Make option in xml to add 3 at once.
                ExtAPILogWriteMessage(str(datetime.now())+" >>> Identifier Slave")
                Identifier[analysisName][elemId]=vals
                # ---- Identifier Slave />
    if not int(result.Analysis.MeshData.ElementById(SelemId).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
        ExtAPILogWriteMessage(str(datetime.now())+" >>> RETURN: Element is not a shell") 
        return []

    if EdgesByElemId[SelemId].Count == 0:
        ExtAPILogWriteMessage(str(datetime.now())+" >>> RETURN: Element is not a weld") 
        return []

    #---- Identifier Master
    maxValue = 0
    for analysisName in ExtAPI.DataModel.AnalysisNames:        
        maxValue = max(maxValue, Identifier[analysisName][SelemId])
    # ExtAPI.Log.WriteMessage("Maxvalue = "+str(maxValue)+"...")
    return maxValue
    #---- Identifier Master />            
    
def selectDisplay(load,prop):
    prop.ClearOptions()
    prop.AddOption("Throat Thickness")
    prop.AddOption("Throat Thickness according to LasE2.mac")
    prop.AddOption(" ")
    prop.AddOption("Scaled Bending over Equivalent (von-Mises) Stress")
    prop.AddOption("Scaled Equivalent (von-Mises) Stress")
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
    prop.AddOption("Envelope")

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
    