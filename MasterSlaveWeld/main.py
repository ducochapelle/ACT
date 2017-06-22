# Philosophy: Calculate everything, design later.
# Compile, Quick, Correct

# T O D O

# 0 Fatigue K factor klopt niet?
# 0.5 make Identifier true int instead of semi string

# BUG: when inserted in an environmen with multiple load cases it will run multiple times. Somehow even if getvalue does nothing, it takes tenths of seconds to run.
# ? Use time to display different results instead of multiple result probes (like weld type or notch case)

# 1  Add 'beta' option, a factor on FW.
# 1.3 enable unselecting elements between analysises. -- zie MultiWeldScaleScope
# -> NOPE, scope all bodies is bugy
# -> also if allowables change, the calculation won't be there
# 1.5 bug: D O E S  N O T  W O R K  I N DesignAssessment
# 1.5 Onendeval data van mesh en tijd en load cases bij houden
# 1.8 clear results moet eigenlijk ook zijn slaves clearen.

# 2 New mesh, new calc, --> wrong database? 
# -> <solver> element?
# -> check d
# 2 bug: Sometimes loadsteps are generated which are not there; probably old tabular data?
# 2 move often used strings to global object for maintainability

# 3 Check if loadcase exists
# 3 possible huge time gain by using ElementValues instead of ElementValue
# 3 import cProfile of profiler, eerst downloaden. Zodat je de echte versie hebt en fatsoenlijk kan profilen.
# 3 fix bug where all the readonly's are on their default value when just inserted; maybe do keypress up and down?

# 4 OOIT automatisch shear factor weld herverdelen

# shear keyoption?  

# System.Diagnostics.Debugger.Break()
if ExtAPI.Context == "Mechanical":
    from datetime import datetime 
    Data = {} 
    WeldElements = [] 
    FaceByElemId = {}
    EdgesByElemId = {}

    # NotchCases = ["W0","W1","W2","K0","K1","K2","K3","K3/4","K4"]
    NotchCases = ["K0","K1","K2","K3"]
    weldThroats = range(3,25,1) + range(25,50,5) + range(50,100,10) + [100,120,150,200]

def ExtAPILogWriteMessage(string):
    # ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> "+str(string))
    # with open("c:\\TEMP\\act.log",'a') as log: # built in logger won't display if there are too many errors or smt.
        # log.write(str(datetime.now())+" >>> "+str(string)+"\n")
    pass
# with open("c:\\TEMP\\act.log", 'w') as log:
    # log.write("If your log was just a text file instead of some 'socket' I could debug the fucking log\n")
    
def DefaultLoadCases():
    LoadCases = {}
    AnalysisNames = ExtAPI.DataModel.AnalysisNames
    for analysisName in AnalysisNames:
        if "Eigenvalue Buckling" in analysisName:
            continue
        try:
            temp = ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.Result("S")
            del temp
        except:
            continue
        LoadCases[analysisName] = {}
        for set in [int(x) for x in ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.ListTimeFreq if x == float(int(x))]:
            LoadCases[analysisName][set]=230
    return LoadCases

def ConvertLoadCases(lc_time):
    """ Steps are used internally, but the user wants to 
    see time. This converts it back or forth.
    Arguments: load cases as a dictionary based on time, analysisName
    Returns: load cases as a dictionary based on load and sub steps
    """
    lc_step = {}
    for analysisName in lc_time:
        d = {}
        ListTimeFreq = ExtAPI.DataModel.AnalysisByName(analysisName).ResultsData.ListTimeFreq
        for time in lc_time[analysisName]:
            d[ListTimeFreq.IndexOf(float(time))+1] = lc_time[analysisName][time]
        lc_step[analysisName]=d
    return lc_step
    
def createMaster(analysis):
    # dit werkt alleen mits je alles in 1 analysis flikkert.
    b = [] # die kut lui geven soms ook een analyse mee die kapot is. Die filter ik er nu dus uit met die try-except.
    for a in ExtAPI.DataModel.AnalysisList:
        try: 
            a.Name
        except:
            continue
        for ro in a.ResultObjects:
            b.append(ro)    
    c = [ro.Properties.GetByName("Identifier").Value for ro in b if ro != None and ro.Name == "Master"]
    for n in [str(n) for n in range(1000)]:
        if n not in c:
            break
    LoadCases = DefaultLoadCases()
    ResultObject = analysis.CreateResultObject("Master")
    for x in ["FatigueCalculation","lc","lcFatigue","CraneGroup","addSlaves","EdgeByElem","EdgeOfModel"]:
        ResultObject.Properties.GetByName(x).ReadOnly = False
    ResultObject.Properties.GetByName("Identifier").Value = n
    ResultObject.Properties.GetByName("lc").Value = str(LoadCases)
    ResultObject.Properties.GetByName("lcFatigue").Value = str(LoadCases)
    ResultObject.Caption = "[{0}] Weld Master".format(n)
    global Data
    Data.pop(n, None) 
    pass
    
def validateLoadCases(result, property):    
    try:
        # als fatigue en dit ding is de fatigue dan moet die maximaal 2 loadcases hebben jeweetoch!
        LoadCases = ConvertLoadCases(eval(property.Value))
        AnalysisNames = [analysisName for analysisName in LoadCases]
        for analysisName in AnalysisNames:
            if not analysisName in ExtAPI.DataModel.AnalysisNames:
                ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" does not exist.")
                property.Value = str(DefaultLoadCases())
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            for resultSet in LoadCases[analysisName]:
                if Analysis.ResultsData.ResultSetCount < resultSet:
                    ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+" has no step "+str(resultSet)+".")
                    property.Value = str(DefaultLoadCases())
                sa = LoadCases[analysisName][resultSet]
                if sa < 10 or sa > 1000:
                    ExtAPI.Application.LogWarning("Analysis "+str(analysisName)+":"+str(resultSet)+" has unexpected allowable: "+str(sa)+" is not between 10 and 1000. Which is already a pretty big range.")
    except:
        ExtAPI.Application.LogWarning("Invalid Load Cases.")
        property.Value = str(DefaultLoadCases())
 
def getNormalDirection(edge,elemId, Mesh):
    # returns the unitvectors alongside the edge, the second one and perpendicular to element
    ExtAPILogWriteMessage("getNormalDirection of edge {0}, element {1}".format(edge.Id, elemId))
    # define two edge nodes (n1, n2)
    edgeNodeIds = set(Mesh.MeshRegionById(edge.Id).NodeIds)
    elemNodeIds = set(Mesh.ElementById(elemId).NodeIds)
    e_e_NodeIds = list(edgeNodeIds & elemNodeIds) # use .pop on set instead ...
    # define vector of edge nodes (A) 
    n1 = Mesh.NodeById(e_e_NodeIds[0])
    n2 = Mesh.NodeById(e_e_NodeIds[1])
    A = ( n2.X - n1.X,
          n2.Y - n1.Y,
          n2.Z - n1.Z)
    A_length = ( A[0]**2+A[1]**2+A[2]**2 ) ** ( 0.5 ) # use sqrt(sum(tuple(x**2 for x in A))) instead...
    A_unit = (A[0]/A_length,
              A[1]/A_length,
              A[2]/A_length)
    # take an additional node which is not in line (n3)
    n3 = Mesh.NodeById((elemNodeIds ^ elemNodeIds & edgeNodeIds).pop())
    # define vector of n1 to n3 (B)
    B = (n3.X - n1.X,
         n3.Y - n1.Y,
         n3.Z - n1.Z)
    # define component of B that is in A (C)
    BdotA_unit = sum([x*y for x,y in zip(B,A_unit)])
    C = [BdotA_unit*A_unit[n] for n in range(3)]
    # define component of B that is not in A, substract C from B (D)
    D = (B[0] - C[0],
         B[1] - C[1],
         B[2] - C[2])
    # unitize it
    D_length = ( D[0]**2+D[1]**2+D[2]**2 ) ** ( 0.5 )
    D_unit = (D[0]/D_length, D[1]/D_length, D[2]/D_length)
    # get third vector perpendicular to the two
    E_unit = (  A_unit[1] * D_unit[2] - A_unit[2] * D_unit[1],
                A_unit[2] * D_unit[0] - A_unit[0] * D_unit[2],
                A_unit[0] * D_unit[1] - A_unit[1] * D_unit[0])
    ExtAPILogWriteMessage("endNormalDirection")
    return (D_unit,A_unit,E_unit)        # normal, length, cross

class Tensor():
    def __init__(self, t, u=None):
        self.Tensor = t
        self.Rotated = None
        self.Normal = None
        self.Long = None
        if u:
            self.Rotate(u)
    def Rotate(self,u):
        self.Rotated = self.mmult(self.mmult(u,self.Tensor),zip(*u))
        self.Normal = self.Rotated[0][0]
        self.Long = self.Rotated[0][1]
    def mmult(self,m1,m2):
        M = [[0,0,0],[0,0,0],[0,0,0]]
        for i in range(len(m2)):
            for j in range(len(m2[0])):
                for k in range(len(m1)):
                    M[i][j] += m1[i][k]*m2[k][j]
                # M[i][j] = round(M[i][j]) # easier to read for debugging
        return M

def createDB(result,step):
    if not 1 == step:
        return
    ExtAPI.Log.WriteMessage("{0} Step: {3}, Start Master Id: {1}, Item: {2}".format(str(datetime.now()), str(result.Properties.GetByName("Identifier").Value), str(result.Caption), str(step)))    
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> CreateDB: Reverse Dictionaries")
    #move all this to a class     
    global WeldElements
    global FaceByElemId                                                                                   
    global EdgesByElemId
    WeldElements = []
    FaceByElemId = {}
    EdgesByElemId = {}

    Mesh = result.Analysis.MeshData
    if WeldElements == []:  # deze loop naar een functie
        ExtAPILogWriteMessage("CreateDB: WeldElementLoop")
        EdgeOfModel = result.Properties.GetByName("EdgeOfModel").Value
        GeoData = result.Analysis.GeoData
        for elem in result.Analysis.MeshData.ElementIds:                                                    
            EdgesByElemId[elem] = []            
        for assembly in result.Analysis.GeoData.Assemblies:
            for part in assembly.Parts:
                for body in part.Bodies:
                    if body.Suppressed == True:
                        continue # if the body is suppressed
                    if not str(body.BodyType) == "GeoBodySheet":
                        continue # if the body is not a shell
                    for face in body.Faces:     #ASLV
                        for edge in face.Edges:     #LSLA
                            if result.Analysis.MeshData.MeshRegionById(edge.Id).NodeCount == 0: 
                                continue # edges which are skipped with meshing
                            if edge.Faces.Count < 2 and EdgeOfModel == "Exclude":    
                                continue # if the edge doesn't have multiple faces; it's not a weld
                            NodeIdsOnEdge = set(Mesh.MeshRegionById(edge.Id).NodeIds)
                            for NodeId in NodeIdsOnEdge:     #NSLL
                                for ConnectedElementId in Mesh.NodeById(NodeId).ConnectedElementIds:    #ESLN 
                                    NodeIdsOfElement = set(Mesh.ElementById(ConnectedElementId).NodeIds)
                                    NodeIdsOfElementOnEdge = NodeIdsOnEdge & NodeIdsOfElement
                                    if len(NodeIdsOfElementOnEdge) == 1 or NodeIdsOfElementOnEdge == NodeIdsOfElement:
                                        continue # if weld direction not clear
                                    if not edge.Id in [x.Id for x in EdgesByElemId[ConnectedElementId]]:
                                        ExtAPILogWriteMessage("CreateDB: EdgesByElemId")
                                        EdgesByElemId[ConnectedElementId].append(edge)
                        for ElementId in Mesh.MeshRegionById(face.Id).ElementIds:
                            ExtAPILogWriteMessage("CreateDB: FaceByElemId")
                            FaceByElemId[ElementId] = face
        ExtAPILogWriteMessage("CreateDB: WeldElements")
        for elem in Mesh.ElementIds:                                       
            if not EdgesByElemId[elem].Count == 0:
                if elem in FaceByElemId:
                    if int(Mesh.ElementById(elem).Type) in (5,6,7,8):   #kTri3, kTri6, kQuad4, kQuad8
                        WeldElements.append(elem)
    Identifier = result.Properties.GetByName("Identifier").Value
    global Data
    Data.pop(Identifier, None)
    if not Identifier in Data:
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> CreateDB: Identifier {0}".format(str(Identifier)))
        Data[Identifier] = {}
        LoadCases = ConvertLoadCases(eval(result.Properties.GetByName("lc").Value))
        EdgeByElem = result.Properties.GetByName("EdgeByElem").Value
        # first halve of calculation is forEachAnalysis(forEachElement) 
        # for performence because .ResultsData takes seconds to load.
        # alternatively, maybe ask all element results in one ElementValues with an s.
        ExtAPILogWriteMessage("CreateDB: Identifier FirstLoop")
        for analysisName in LoadCases: 
            Data[Identifier][analysisName] = {}
            Analysis = ExtAPI.DataModel.AnalysisByName(analysisName)
            for resultSet in LoadCases[analysisName]:
                Data[Identifier][analysisName][resultSet] = {}
                #maybe ask stress of the nodes who are actual on the line?? 
                #I think the element Ids can be easily queried...
                Analysis.ResultsData.CurrentResultSet = resultSet
                resultStress = Analysis.ResultsData.Result("S")
                allowableStress = LoadCases[analysisName][resultSet]
                
                for elemId in WeldElements:   
                    ExtAPILogWriteMessage("CreateDB: Analysis:{0}{1}, elemId:{2}".format(analysisName, resultSet, str(elemId)))
                    Data[Identifier][analysisName][resultSet][elemId]={"None":{"Static":{}}
                                                                      ,"2xFW":{"Static":{}}
                                                                      ,"2xPP":{"Static":{}}
                                                                      ,"1xPP":{"Static":{}}}
                    DeepData = Data[Identifier][analysisName][resultSet][elemId]["None"]["Static"]
                    plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
                    # stresses are no longer abs() -> check rest of this indent careful for consequences
                    sxxTotal = resultStress.ElementValue(elemId,"X")
                    syyTotal = resultStress.ElementValue(elemId,"Y")
                    szzTotal = resultStress.ElementValue(elemId,"Z")
                    sxyTotal = resultStress.ElementValue(elemId,"XY")
                    syzTotal = resultStress.ElementValue(elemId,"YZ")
                    sxzTotal = resultStress.ElementValue(elemId,"XZ")
                    
                    #0:4 = top plane, 4:8 = bottom plane, 8:12 = mid plane
                    sxx = sum(sxxTotal[0:4],0.0)/4  
                    syy = sum(syyTotal[0:4],0.0)/4  
                    szz = sum(szzTotal[0:4],0.0)/4  
                    sxy = sum(sxyTotal[0:4],0.0)/4  
                    syz = sum(syzTotal[0:4],0.0)/4  
                    sxz = sum(sxzTotal[0:4],0.0)/4  
                    TensorTop = Tensor( ((  sxx,    sxy,    sxz)
                                        ,(  sxy,    syy,    syz)
                                        ,(  sxz,    syz,    szz))  )

                    sxx = sum(sxxTotal[8:12],0.0)/4
                    syy = sum(syyTotal[8:12],0.0)/4
                    szz = sum(szzTotal[8:12],0.0)/4
                    sxy = sum(sxyTotal[8:12],0.0)/4
                    syz = sum(syzTotal[8:12],0.0)/4
                    sxz = sum(sxzTotal[8:12],0.0)/4
                    TensorMid = Tensor( ((  sxx,    sxy,    sxz)
                                        ,(  sxy,    syy,    syz)
                                        ,(  sxz,    syz,    szz))  )
                    ExtAPILogWriteMessage("CreateDB: Got Stresses")
                    
                    # Determine the edge to calculate with; for elements can be with multiple edges
                    # Contains duplicate calculations to be sure the worst-case edge is taken into account
                    if EdgesByElemId[elemId].Count == 1:
                        Edge = EdgesByElemId[elemId][0]
                    else:
                        ExtAPILogWriteMessage("CreateDB: EnvelopeStart")    # do this more thoroughyfully (including moment). OHja? EHEH PP,FW wie ruled dan?
                        if EdgeByElem == "Envelope":
                            CombinedSqStress = 0
                            for edge in EdgesByElemId[elemId]:
                                units = getNormalDirection(edge, elemId, Mesh) 
                                TensorMid.Rotate(units)
                                NormalStress    = TensorMid.Normal
                                LongStress      = TensorMid.Long
                                newCombinedSqStress = NormalStress**2 + 3*LongStress**2
                                if newCombinedSqStress >= CombinedSqStress:
                                    EdgeWithMaxWeldThickness = edge
                                    CombinedSqStress = newCombinedSqStress
                            Edge = EdgeWithMaxWeldThickness
                        elif EdgeByElem == "High Edge Ids":
                            Edge = max([x.Id for x in EdgesByElemId[elemId]])
                        elif EdgeByElem == "Low Edge Ids":
                            Edge = min([x.Id for x in EdgesByElemId[elemId]])
                        elif EdgeByElem == "Random Edge Id":
                            Edge = EdgesByElemId[elemId][0]
                        ExtAPILogWriteMessage("CreateDB: EnvelopeStop")
                            
                    # Calculate the actual weld
                    units = getNormalDirection(Edge, elemId, Mesh)
                    # don't use edge but edge nodes? Check if edge is GeoCurveLine or GeoCurveCircle
                    TensorMid.Rotate(units)
                    TensorTop.Rotate(units)
                    
                    NormalStress    = TensorMid.Normal
                    LongStress      = TensorMid.Long
                    EquivalentStress = ( NormalStress**2+3*LongStress**2 )**(0.5) 
                    weldThickness = ( EquivalentStress * plateThickness ) / (allowableStress)
                    DeepData["Equivalent (von-Mises) Stress"]           = EquivalentStress
                    DeepData["Axial Stress"]                            = NormalStress #abs?
                    DeepData["Longitudinal Stress"]                     = LongStress #abs?
                    DeepData["Throat Thickness according to LasE2.mac"] = weldThickness
                    
                    NormalForce = NormalStress * plateThickness
                    LongForce = LongStress * plateThickness
                    Moment = (TensorTop.Normal - TensorMid.Normal) * plateThickness**2 * 6**(-1) 
                    DeepData["Axial Force mm^-1"]           = NormalForce
                    DeepData["Longitudinal Moment mm^-1"]   = Moment
                    DeepData["Longitudinal Force mm^-1"]    = LongForce
                    ExtAPILogWriteMessage("CreateDB: Calculated Forces")
                    # Get unit vectors for debuggin purpose
                    unitNormalScalar = round(abs(9*units[0][0]))*100 \
                                     + round(abs(9*units[0][1]))*10 \
                                     + round(abs(9*units[0][2]))
                    unitCrossScalar  = round(abs(9*units[2][0]))*100 \
                                     + round(abs(9*units[2][1]))*10 \
                                     + round(abs(9*units[2][2]))
                    unitLongScalar   = round(abs(9*units[1][0]))*100 \
                                     + round(abs(9*units[1][1]))*10 \
                                     + round(abs(9*units[1][2]))
                    DeepData["Axial Unit Vector"]           = unitNormalScalar
                    DeepData["Cross Unit Vector"]           = unitCrossScalar
                    DeepData["Longitudinal Unit Vector"]    = unitLongScalar
                    for weldType in ["2xFW","2xPP","1xPP"]:
                        allowableStress = LoadCases[analysisName][resultSet]
                        notchCase = "Static"
                        for a in weldThroats:       
                            if weldType == "2xFW":
                                sigma_n = tau_n = NormalForce * 1.414 / (4 * a)
                                tau_l           = LongForce / (2 * a)
                                b               = plateThickness + 0.667 * a * 1.414
                                sigma_m         = tau_m = Moment / (1.414 * a * b)
                            elif weldType == "2xPP":
                                sigma_n         = NormalForce / (2 * a)
                                tau_l           = LongForce / (2 * a)
                                b               = plateThickness - a
                                if b == 0: continue
                                sigma_m = Moment / (a * b)  
                                tau_n = tau_m = 0
                            elif weldType == "1xPP":
                                sigma_n         = NormalForce / a
                                tau_l           = LongForce / a
                                sigma_m         = Moment / (0.250 * a**2)
                                tau_n = tau_m   = 0
                            sigma_vm1   = sqrt((sigma_n + sigma_m)**2 + 3*(tau_l)**2 + 3*(tau_n + tau_m)**2)
                            sigma_vm2   = sqrt((sigma_n - sigma_m)**2 + 3*(tau_l)**2 + 3*(tau_n - tau_m)**2)
                            sigma_vm    = max(sigma_vm1, sigma_vm2)
                            if sigma_vm < allowableStress: 
                                axialStressInWeld   = sigma_n + tau_n   # similar to sqrt(sigma_n**2 + 3*tau_n**2) given that sigma == tau or 0 == tau
                                bendingStressInWeld = sigma_m + tau_m   # similar to sqrt(sigma_m**2 + 3*tau_m**2) given that sigma == tau or 0 == tau
                                longStressInWeld    = 1.732*tau_l       # similar to sqrt(             3*tau_l**2)
                                DeepData = Data[Identifier][analysisName][resultSet][elemId][weldType][notchCase]
                                DeepData["Throat Thickness"]                       = a
                                DeepData["Scaled Axial Stress"]                    = axialStressInWeld  
                                DeepData["Scaled Bending Stress"]                  = bendingStressInWeld
                                DeepData["Scaled Longitudinal Stress"]             = longStressInWeld   
                                DeepData["Scaled Equivalent (von-Mises) Stress"]   = sigma_vm
                                if sigma_vm < 1:
                                    sigma_vm = 1
                                DeepData["Scaled Axial over Equivalent (von-Mises) Stress"]        = axialStressInWeld / sigma_vm
                                DeepData["Scaled Bending over Equivalent (von-Mises) Stress"]      = bendingStressInWeld / sigma_vm
                                DeepData["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = longStressInWeld / sigma_vm
                                ExtAPILogWriteMessage("CreateDB: Found a < 200")
                                break
                        else:
                            DeepData = Data[Identifier][analysisName][resultSet][elemId][weldType][notchCase]
                            a = 999
                            DeepData["Throat Thickness"]                       = a
                            DeepData["Scaled Axial Stress"]                    = a
                            DeepData["Scaled Bending Stress"]                  = a
                            DeepData["Scaled Longitudinal Stress"]             = a
                            DeepData["Scaled Equivalent (von-Mises) Stress"]   = a
                            DeepData["Scaled Axial over Equivalent (von-Mises) Stress"]        = a
                            DeepData["Scaled Bending over Equivalent (von-Mises) Stress"]      = a
                            DeepData["Scaled Longitudinal over Equivalent (von-Mises) Stress"] = a
                            ExtAPILogWriteMessage("CreateDB: Found a > 200")
                            # ExtAPILogWriteMessage("NormalStress {} LongStress {} TensorTop.Normal {} TensorMid.Normal {} plateThickness {}".format(NormalStress,LongStress,TensorTop.Normal,TensorMid.Normal,plateThickness)) # FFFUUUUU
    ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> CreateDB: Finished Static Without Errors")
    if not result.Properties.GetByName("FatigueCalculation").Value == "No": 
        #static calculation performence should not be affected at any rate by fatigue calculation
        CalculateFatigue(result.Properties, Data[Identifier])
        ExtAPI.Log.WriteMessage(str(datetime.now())+" >>> CreateDB: Finished Fatigue Without Errors")

def CalculateFatigue(Properties, DataId):
    LoadCases = ConvertLoadCases(eval(Properties.GetByName("lcFatigue").Value))
    CraneGroup = Properties.GetByName("CraneGroup").Value
    for elemId in WeldElements:
        plateThickness = FaceByElemId[elemId].Body.Thickness * 1000   # to mm, this is so ugly
        for weldType in ["2xFW","2xPP","1xPP"]:
            a = 0
            for analysisName in LoadCases:  # get maximum throat thickness from static calculations
                for resultSet in LoadCases[analysisName]:
                    a = max(a,DataId[analysisName][resultSet][elemId][weldType]["Static"]["Throat Thickness"])
            WeldStressDirections = []
            for analysisName in LoadCases:  # get stress directions
                for resultSet in LoadCases[analysisName]:
                    DeepData = DataId[analysisName][resultSet][elemId]["None"]["Static"]
                    NormalForce = DeepData["Axial Force mm^-1"]
                    LongForce   = DeepData["Longitudinal Force mm^-1"]
                    Moment      = DeepData["Longitudinal Moment mm^-1"]
                    ExtAPILogWriteMessage("Analysis: {0}, resultSet: {1}, WeldType: {2}, Normal: {3}, Long: {4}, Moment: {5}".format(analysisName, resultSet, weldType, NormalForce, LongForce, Moment))
                    if weldType == "2xFW":
                        sigma_n = tau_n = NormalForce * 1.414 / (4 * a)
                        tau_l           = LongForce / (2 * a)
                        b               = plateThickness + 0.667 * a * 1.414
                        sigma_m         = tau_m = Moment / (1.414 * a * b)
                    elif weldType == "2xPP":
                        sigma_n         = NormalForce / (2 * a)
                        tau_l           = LongForce / (2 * a)
                        b               = plateThickness - a
                        if b == 0: continue
                        sigma_m = Moment / (a * b)  
                        tau_n = tau_m = 0
                    elif weldType == "1xPP":
                        sigma_n         = NormalForce / a
                        tau_l           = LongForce / a
                        sigma_m         = Moment / (0.250 * a**2)
                        tau_n = tau_m   = 0
                    axialStressInWeld   = sigma_n + tau_n # similar to sqrt(sigma_n**2 + 3*tau_n**2) given that sigma == tau or 0 == tau
                    bendingStressInWeld = sigma_m + tau_m # similar to sqrt(sigma_m**2 + 3*tau_m**2) given that sigma == tau or 0 == tau
                    longStressInWeld    = 1.732*tau_l # similar to sqrt(             3*tau_l**2)
                    # get both (ax - m, ax + m, long), define K factor
                    WeldStressDirections.append({"ax1":axialStressInWeld - bendingStressInWeld, \
                                                 "ax2":axialStressInWeld + bendingStressInWeld, \
                                                 "lng":longStressInWeld}) 
            ExtAPILogWriteMessage(WeldStressDirections)
            K = CalculateK(WeldStressDirections)   
            for analysisName in LoadCases: 
                for resultSet in LoadCases[analysisName]: 
                    # for future reference. Dit is wel een beetje erin 
                    # gepropt hoor.. zodat er maar geen code aangepast 
                    # hoeft te worden. Dit is niet echt het meest voor 
                    # de hand liggende plekje van iets. Gewoon niet logisch
                    DataId[analysisName][resultSet][elemId][weldType]["Static"]["K-factor"]=K 
                    for notchCase in NotchCases:
                        allowableStress = FatigueStresses.Get(CraneGroup,str(K),notchCase)
                        DeepData = DataId[analysisName][resultSet][elemId]["None"]["Static"]
                        NormalForce = DeepData["Axial Force mm^-1"]
                        LongForce   = DeepData["Longitudinal Force mm^-1"]
                        Moment      = DeepData["Longitudinal Moment mm^-1"]
                        for a in weldThroats:       
                            if weldType == "2xFW":
                                sigma_n = tau_n = NormalForce * 1.414 / (4 * a)
                                tau_l           = LongForce / (2 * a)
                                b               = plateThickness + 0.670 * a * 1.414
                                sigma_m         = tau_m = Moment / (1.414 * a * b)
                            elif weldType == "2xPP":
                                sigma_n         = NormalForce / (2 * a)
                                tau_l           = LongForce / (2 * a)
                                b               = plateThickness - a
                                if b == 0: continue
                                sigma_m = Moment / (a * b)  
                                tau_n = tau_m = 0
                            elif weldType == "1xPP":
                                sigma_n         = NormalForce / a
                                tau_l           = LongForce / a
                                sigma_m         = Moment / (0.250 * a**2)
                                tau_n = tau_m   = 0
                            sigma_vm1   = sqrt((sigma_n + sigma_m)**2 + 3*(tau_l)**2 + 3*(tau_n + tau_m)**2)
                            sigma_vm2   = sqrt((sigma_n - sigma_m)**2 + 3*(tau_l)**2 + 3*(tau_n - tau_m)**2)
                            sigma_vm    = max(sigma_vm1, sigma_vm2)
                            if sigma_vm < allowableStress: 
                                DataId[analysisName][resultSet][elemId][weldType][notchCase]={"Throat Thickness":a}
                                ExtAPILogWriteMessage("CreateDB: Found a < 200 Fatigue")
                                break
                        else:
                            DataId[analysisName][resultSet][elemId][weldType][notchCase]={"Throat Thickness":999}
                            ExtAPILogWriteMessage("CreateDB: Found a > 200 Fatigue")
                        
def CalculateK(WeldStressDirections):   
    """ Calulates K from three stress directions in two load cases.
    
    >>> CalculateK([{"ax1":10.0,"ax2":10.0,"lng":10.0},{"ax1":10.0,"ax2":10.0,"lng":10.0}])
    1.0
    >>> CalculateK([{"ax1":0.0,"ax2":0.0,"lng":0.0},{"ax1":10.0,"ax2":10.0,"lng":10.0}])
    0.0
    >>> CalculateK([{"ax1":-10.0,"ax2":-10.0,"lng":-10.0},{"ax1":10.0,"ax2":10.0,"lng":10.0}])
    -1.0
    >>> CalculateK([{"ax1":10.0,"ax2":10.0,"lng":1.0},{"ax1":10.0,"ax2":10.0,"lng":-1.0}])
    0.8
    >>> CalculateK([{"ax1":-10.0,"ax2":-10.0,"lng":10.0},{"ax1":10.0,"ax2":10.0,"lng":10.0}])
    0.0
    >>> CalculateK([{"ax1":-10.0,"ax2":-10.0,"lng":0.0},{"ax1":10.0,"ax2":10.0,"lng":0.0}])
    -1.0
    >>> CalculateK([{"ax1":0.0,"ax2":0.0,"lng":0.0},{"ax1":0.0,"ax2":0.0,"lng":0.0}])
    1.0
    >>> CalculateK([{'lng': 74.3234926291666, 'ax1': -109.87699064127602, 'ax2': 180.72234464518229}, {'lng': 16.771171839519191, 'ax1': -71.454549560546866, 'ax2': 218.79426879882811}])
    0.5
    """
    # todo: check lcFatigue 
    # todo: en als het allmeaal nul is? even goed doordenken! EN dus niet gewoon alles 1 maken en ervanuitgaan dat dat dan wel goed zit
    ax1s = [d["ax1"] for d in WeldStressDirections]
    ax2s = [d["ax2"] for d in WeldStressDirections]
    lngs = [d["lng"] for d in WeldStressDirections]
    maxAx1 = max([abs(x) for x in ax1s])
    maxAx2 = max([abs(x) for x in ax2s])
    maxLng = max([abs(x) for x in lngs])
    AxFactor1 = maxAx1/(maxAx1+maxLng) if maxAx1+maxLng != 0 else 1
    AxFactor2 = maxAx2/(maxAx2+maxLng) if maxAx2+maxLng != 0 else 1
    K_ax1 = maxabs_over_minabs(ax1s)
    K_ax2 = maxabs_over_minabs(ax2s)
    K_lng = maxabs_over_minabs(lngs)
    K_1 = K_ax1 * AxFactor1 + K_lng * (1 - AxFactor1)  
    K_2 = K_ax2 * AxFactor2 + K_lng * (1 - AxFactor2)  
    K = round(min(K_1, K_2),1)
    K = 0.0 if K == -0.0 else K # ...
    # print str([str(x) for x in [K, ax1s, ax2s, lngs, maxAx1, maxAx2, maxLng, AxFactor1, AxFactor2, K_ax1, K_ax2, K_lng, K_1, K_2]])
    return K
    
def maxabs_over_minabs(x):
    """ Returns the fraction of the absolutely smallest over largest value.
    
    >>> maxabs_over_minabs((1.0,2.0))
    0.5
    >>> maxabs_over_minabs((-1.0,2.0))
    -0.5
    >>> maxabs_over_minabs((1.0,-2.0))
    -0.5
    >>> maxabs_over_minabs((-1.0,-2.0))
    0.5
    >>> maxabs_over_minabs((-0.0,-2.0))
    0.0
    >>> maxabs_over_minabs((-2.0,-2.0))
    1.0
    >>> maxabs_over_minabs((-2.0,2.0))
    -1.0
    >>> maxabs_over_minabs((-2.0,0.4))
    -0.2
    """
    return min(x)/max(x) if abs(min(x)) < max(x) and max(x) != 0 else max(x)/min(x) if min(x) != 0 else 1
    
def getValue(result,elemId):
    if sp == None:
        return [0.0] # dit kost die kut al tienden van secondes. 
    if elemId not in WeldElements:
        return []
    Identifier      = sp.Identifier      
    Type            = sp.Type            
    Item            = sp.Item            
    WeldType        = sp.WeldType        
    NotchCase       = sp.NotchCase       
    StaticLoadCases = sp.StaticLoadCases 
    CraneGroup      = sp.CraneGroup      
    FatigueLoadCases= sp.FatigueLoadCases

    if not Identifier in Data:
        return [0.0]

    # Ok dit is een ingewikkeld stuk... maar het kost echt nul computation time... superraar.
    vals = []
    sals = []
    if Item == "Scaled Longitudinal over Equivalent (von-Mises) Stress":
        for analysisName in StaticLoadCases:
            for resultSet in StaticLoadCases[analysisName]:
                Deep = Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase]
                vals.append( Deep["Throat Thickness"])                      # regular a
                sals.append((Deep["Throat Thickness"] * (1 - Deep[Item])))  # part of the a which cannot be averaged
        return [ abs( 1 - max(sals) / max(vals) ) ]
    elif Type == "No":
        for analysisName in StaticLoadCases:
            for resultSet in StaticLoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
    elif Item == "K-factor":
        for analysisName in FatigueLoadCases:
            for resultSet in FatigueLoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
    elif Type == "Fat-o-Stat" and NotchCase == "Static":
        for analysisName in StaticLoadCases:
            for resultSet in StaticLoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
    elif Type == "Combined":
        for analysisName in StaticLoadCases:
            for resultSet in StaticLoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType]["Static"][Item])
                if analysisName in FatigueLoadCases and resultSet in FatigueLoadCases[analysisName]:
                    vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
    elif Type == "Seperate":
        LoadCases = StaticLoadCases if NotchCase == "Static" else FatigueLoadCases
        for analysisName in LoadCases:
            for resultSet in LoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
    elif Type == "Fat-o-Stat":
        for analysisName in FatigueLoadCases:
            for resultSet in FatigueLoadCases[analysisName]:
                vals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType][NotchCase][Item])
                sals.append(Data[Identifier][analysisName][resultSet][elemId][WeldType]["Static"][Item])
        if max(vals) == 5.0: 
            return [0.0] # te kort door de bocht?
        else:
            return [float(max(vals)) / float(max(sals))]
    return [max(vals)*1.0]

def selectEdgeByElem(load,prop):
    prop.ClearOptions()
    prop.AddOption("Envelope")
    prop.AddOption("High Edge Ids")
    prop.AddOption("Low Edge Ids")
    prop.AddOption("Random Edge Id")
    
def ClearResults(result, prop):
    ExtAPI.Log.WriteMessage("Clear Results...")
    Identifier = result.Properties.GetByName("Identifier").Value
    for ResultObject in [ro for ro in result.Analysis.ResultObjects if ro.Name == "Slave" and ro.Properties.GetByName("Identifier").Value == Identifier]:
        pass # mughaga just delete and recreate the items? no fakckkk
    global WeldElements
    global FaceByElemId
    global EdgesByElemId
    WeldElements = []
    FaceByElemId = {}
    EdgesByElemId = {}
    global Data
    Data.pop(Identifier, None)
    
def IsValid(foo,bar):
    return True

def selectFatigueCalculation(result,prop):
    prop.ClearOptions()
    prop.AddOption("No")
    prop.AddOption("Combined")
    prop.AddOption("Seperate")
    prop.AddOption("Fat-o-Stat")
    
def validateFatigueCalculation(result,prop):    
    if prop.Value == "No":
        result.Properties.GetByName("lcFatigue").Visible = False
        result.Properties.GetByName("CraneGroup").Visible = False
    else:
        result.Properties.GetByName("lcFatigue").Visible = True
        result.Properties.GetByName("CraneGroup").Visible = True
    
def selectEdgeOfModel(result,prop):
    prop.ClearOptions()
    prop.AddOption("Include")
    prop.AddOption("Exclude")
    
def selectCraneGroup(result,prop):
    prop.ClearOptions()
    for o in ["B"+str(n) for n in range(7)]: 
        prop.AddOption(o)
    
def validateCraneGroup(result,prop):
    pass
    
def addFatigueSlaves(result, prop): 
    # ipv al die shit mee geven, gewoon checken bij de master?
    # De K factor ook hier plotten
    Type = result.Properties.GetByName("FatigueCalculation").Value
    global NotchCases
    MProp = result.Properties
    for wt in ["2xFW", "2xPP", "1xPP"]:
        for nc in ["Static"] + NotchCases:
            ResultObject = result.Analysis.CreateResultObject("Slave")
            SProp = ResultObject.Properties
            if nc == "Static": SProp.GetByName("Item").Value  = "K-factor"
            SProp.GetByName("Geometry").Value               = MProp.GetByName("Geometry").Value
            SProp.GetByName("Identifier").Value             = MProp.GetByName("Identifier").Value
            SProp.GetByName("FatigueCalculation").Value     = MProp.GetByName("FatigueCalculation").Value 
            SProp.GetByName("lc").Value                     = str(ConvertLoadCases(eval(MProp.GetByName("lc").Value)))
            SProp.GetByName("lcFatigue").Value              = MProp.GetByName("lcFatigue").Value
            SProp.GetByName("CraneGroup").Value             = MProp.GetByName("CraneGroup").Value
            SProp.GetByName("NotchCase").Value              = nc
            SProp.GetByName("WeldType").Value               = wt
            Caption = ""
            Caption += "["+SProp.GetByName("Identifier").Value+"]" + " "
            Caption += SProp.GetByName("WeldType").Value + " "
            if nc == "Static":
                Caption += "K-factor"
            else:
                Caption += SProp.GetByName("NotchCase").Value
            ResultObject.Caption = Caption
            
def addStaticSlaves(result, prop): 
    Type = result.Properties.GetByName("FatigueCalculation").Value
    MProp = result.Properties
    for i in range(2):
        for wt in ["2xFW", "2xPP", "1xPP"]:
            ResultObject = result.Analysis.CreateResultObject("Slave")
            SProp = ResultObject.Properties
            SProp.GetByName("Geometry").Value               = MProp.GetByName("Geometry").Value
            SProp.GetByName("Identifier").Value             = MProp.GetByName("Identifier").Value
            SProp.GetByName("FatigueCalculation").Value     = MProp.GetByName("FatigueCalculation").Value 
            if i == 1: # needs range(2)...
                SProp.GetByName("Item").Value               = "Scaled Longitudinal over Equivalent (von-Mises) Stress"
            SProp.GetByName("lc").Value                     = str(ConvertLoadCases(eval(MProp.GetByName("lc").Value)))
            SProp.GetByName("lcFatigue").Visible            = False
            SProp.GetByName("CraneGroup").Visible           = False
            SProp.GetByName("NotchCase").Value              = "Static"
            SProp.GetByName("WeldType").Value               = wt
            Caption = ""
            Caption += "["+SProp.GetByName("Identifier").Value+"]" + " "
            Caption += SProp.GetByName("WeldType").Value + " "
            Caption += SProp.GetByName("NotchCase").Value
            if i == 1:
                Caption += " Shear Factor"
            ResultObject.Caption = Caption

def addSlaves(result,prop):
    Type = result.Properties.GetByName("FatigueCalculation").Value
    if Type != "No" and result.Properties.GetByName("lcFatigue").Value == "{}":
        ExtAPI.Application.LogWarning("Fatigue load cases empty")
        return
    if Type == "No":
        addStaticSlaves(result,prop)
    elif Type == "Combined":
        addFatigueSlaves(result,prop)
    elif Type == "Seperate":
        addStaticSlaves(result,prop)
        addFatigueSlaves(result,prop)
    elif Type == "Fat-o-Stat":
        addStaticSlaves(result,prop)
        addFatigueSlaves(result,prop)
    for x in ["FatigueCalculation","lc","lcFatigue","CraneGroup","addSlaves","EdgeByElem","EdgeOfModel"]:
        result.Properties.GetByName(x).ReadOnly = True
    result.Caption = result.Properties.GetByName("Identifier").Value + " " + result.Properties.GetByName("FatigueCalculation").Value # DOES NOT COMPUTE
    
def addMoreSlaves(result,prop):      # this option is intended for debuggin only, only use this when you read the code!
    a = ["Throat Thickness", \
        "Scaled Axial Stress", \
        "Scaled Bending Stress", \
        "Scaled Longitudinal Stress", \
        "Scaled Equivalent (von-Mises) Stress", \
        "Scaled Longitudinal over Equivalent (von-Mises) Stress"]
    b = ["Equivalent (von-Mises) Stress", \
        "Axial Stress", \
        "Longitudinal Stress", \
        "Throat Thickness according to LasE2.mac", \
        "Axial Force mm^-1", \
        "Longitudinal Moment mm^-1", \
        "Longitudinal Force mm^-1", \
        "Axial Unit Vector", \
        "Cross Unit Vector", \
        "Longitudinal Unit Vector"]        
    MProp = result.Properties    
    for wt in ["None","2xFW","2xPP","1xPP"]:
        if wt == "None":
            items = b
        else:
            items = a
        for item in items:
            ResultObject = result.Analysis.CreateResultObject("Slave")
            SProp = ResultObject.Properties
            SProp.GetByName("Item").Value                   = item
            Sprop.GetByName("Geometry").Value               = MProp.GetByName("Geometry").Value
            SProp.GetByName("Identifier").Value             = MProp.GetByName("Identifier").Value
            SProp.GetByName("FatigueCalculation").Value     = MProp.GetByName("FatigueCalculation").Value 
            SProp.GetByName("lc").Value                     = str(ConvertLoadCases(eval(MProp.GetByName("lc").Value)))
            SProp.GetByName("lcFatigue").Visible            = False
            SProp.GetByName("CraneGroup").Visible           = False
            SProp.GetByName("NotchCase").Value              = "Static"
            SProp.GetByName("WeldType").Value               = wt
            Caption = ""
            Caption += "["+SProp.GetByName("Identifier").Value+"]" + " "
            Caption += SProp.GetByName("WeldType").Value + " "
            Caption += SProp.GetByName("NotchCase").Value + " "
            Caption += SProp.GetByName("Item").Value
            ResultObject.Caption = Caption
    
class cFatigueStresses:
    def __init__(self):
        # tab seperated, 8 spaces indented
        # but why is a K1 better than a K0 when K is high?
        self.String = """
        B0
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    669    432    367    511    511    429    307    245    184
        -0.9    697    450    382    532    532    447    319    256    192
        -0.8    727    469    399    555    555    467    333    267    200
        -0.7    760    491    417    581    581    488    348    279    209
        -0.6    797    514    437    608    608    511    365    292    219
        -0.5    836    540    459    639    639    537    383    307    230
        -0.4    880    568    483    672    672    565    403    323    242
        -0.3    929    600    510    710    710    596    426    341    255
        -0.2    984    635    540    751    751    631    451    361    270
        -0.1    999    675    573    798    798    671    479    383    287
        0.0    999    719    611    852    852    715    511    409    307
        0.1    999    785    667    912    929    780    557    446    334
        0.2    999    863    734    983    999    858    613    491    368
        0.3    999    959    815    999    999    954    681    545    409
        0.4    999    999    917    999    999    999    767    613    460
        0.5    999    999    999    999    999    999    876    701    525
        0.6    999    999    999    999    999    999    999    818    613
        0.7    999    999    999    999    999    999    999    981    736
        0.8    999    999    999    999    999    999    999    999    919
        0.9    999    999    999    999    999    999    999    999    999
        1.0    999    999    999    999    999    999    999    999    999
        B1
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    498    321    273    371    371    312    223    178    134
        -9.0    519    335    284    387    387    325    232    186    139
        -8.0    541    349    297    403    403    339    242    194    145
        -7.0    566    365    310    422    422    354    253    202    152
        -6.0    593    383    325    442    442    371    265    212    159
        -5.0    623    402    341    464    464    390    278    223    167
        -4.0    655    423    359    488    488    410    293    234    176
        -3.0    692    446    379    515    515    433    309    247    186
        -2.0    733    473    402    546    546    458    327    262    196
        -1.0    778    502    427    580    580    487    348    278    209
        0.0    830    536    455    619    619    519    371    297    223
        0.1    889    584    497    663    675    567    405    324    243
        0.2    958    643    546    714    742    623    445    356    267
        0.3    999    714    607    773    825    693    495    396    297
        0.4    999    803    683    843    928    779    557    445    334
        0.5    999    918    780    928    999    891    636    509    382
        0.6    999    999    910    999    999    999    742    594    445
        0.7    999    999    999    999    999    999    891    712    534
        0.8    999    999    999    999    999    999    999    890    668
        0.9    999    999    999    999    999    999    999    999    999
        1.0    999    999    999    999    999    999    999    999    999
        B2
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    371    239    203    278    270    226    162    129    97
        -9.0    386    249    212    289    281    236    168    135    101
        -8.0    403    260    221    302    293    246    176    141    105
        -7.0    421    272    231    316    306    257    184    147    110
        -6.0    441    285    242    331    321    270    192    154    115
        -5.0    463    299    254    347    337    283    202    162    121
        -4.0    488    315    268    366    355    298    213    170    128
        -3.0    515    332    282    386    374    314    225    180    135
        -2.0    545    352    299    409    396    333    238    190    143
        -1.0    579    374    318    434    421    354    253    202    152
        0.0    618    399    339    463    449    377    269    216    162
        0.1    662    435    370    496    490    412    294    235    176
        0.2    713    478    407    534    539    453    323    259    194
        0.3    772    532    452    579    599    503    359    288    216
        0.4    842    598    508    632    674    566    404    324    242
        0.5    927    683    581    695    770    647    462    370    277
        0.6    999    797    678    772    898    755    539    431    323
        0.7    999    957    813    868    999    906    647    518    388
        0.8    999    999    999    999    999    999    808    647    485
        0.9    999    999    999    999    999    999    999    999    970
        1.0    999    999    999    999    999    999    999    999    999
        B3
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    276    178    151    234    196    164    118    94    71
        -0.9    288    186    158    244    204    171    122    98    73
        -0.8    300    194    164    254    213    179    128    102    77
        -0.7    314    202    172    266    223    187    134    107    80
        -0.6    329    212    180    279    233    196    140    112    84
        -0.5    345    223    189    293    245    206    147    118    88
        -0.4    363    234    199    308    258    216    155    124    93
        -0.3    383    247    210    325    272    228    163    131    98
        -0.2    406    262    223    344    288    242    173    138    104
        -0.1    431    278    236    366    306    257    184    147    110
        0.0    460    297    252    390    326    274    196    157    117
        0.1    493    324    275    418    356    299    214    171    128
        0.2    531    356    303    450    392    329    235    188    141
        0.3    575    396    336    487    435    365    261    209    157
        0.4    627    445    378    532    490    411    294    235    176
        0.5    690    509    432    585    559    470    336    269    201
        0.6    767    594    504    650    653    548    392    313    235
        0.7    862    712    605    731    783    658    470    376    282
        0.8    986    890    756    836    979    822    587    470    352
        0.9    999    999    999    999    999    999    999    940    705
        1.0    999    999    999    999    999    999    999    999    999
        B4
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    205    133    113    197    142    119    85    68    51
        -0.9    214    138    117    205    148    124    89    71    53
        -0.8    223    144    123    214    155    130    93    74    56
        -0.7    233    151    128    224    162    136    97    78    58
        -0.6    245    158    134    235    169    142    102    81    61
        -0.5    257    166    141    246    178    149    107    85    64
        -0.4    270    174    148    259    187    157    112    90    67
        -0.3    285    184    157    274    197    166    118    95    71
        -0.2    302    195    166    290    209    176    125    100    75
        -0.1    321    207    176    308    222    187    133    107    80
        0.0    342    221    188    328    237    199    142    114    85
        0.1    367    241    205    352    259    217    155    124    93
        0.2    395    265    225    379    284    239    171    136    102
        0.3    428    294    250    411    316    265    190    152    114
        0.4    467    331    282    448    355    299    213    171    128
        0.5    514    379    322    493    406    341    244    195    146
        0.6    571    442    376    547    474    398    284    227    171
        0.7    642    530    451    616    569    478    341    273    205
        0.8    734    662    563    704    711    597    426    341    256
        0.9    999    999    999    999    999    999    853    682    512
        1.0    999    999    999    999    999    999    999    999    999
        B5
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    155    100    85    155    103    87    62    50    37
        -0.9    161    104    89    161    108    90    65    52    39
        -0.8    168    109    92    168    112    94    67    54    40
        -0.7    176    114    97    176    117    99    70    56    42
        -0.6    185    119    101    185    123    103    74    59    44
        -0.5    194    125    106    194    129    108    78    62    47
        -0.4    204    132    112    204    136    114    82    65    49
        -0.3    215    139    118    215    143    120    86    69    52
        -0.2    228    147    125    228    152    128    91    73    55
        -0.1    242    156    133    242    161    135    97    78    58
        0.0    258    167    142    258    172    145    103    83    62
        0.1    277    182    155    277    188    158    113    90    68
        0.2    298    200    170    298    207    173    124    99    74
        0.3    323    222    189    323    230    193    138    110    83
        0.4    352    250    212    352    258    217    155    124    93
        0.5    387    286    243    387    295    248    177    142    106
        0.6    431    333    283    431    344    289    207    165    124
        0.7    484    400    340    484    413    347    248    198    149
        0.8    554    500    425    554    516    433    310    248    186
        0.9    999    999    850    999    999    867    620    496    372
        1.0    999    999    999    999    999    999    999    999    999
        B6
        K    W0    W1    W2    K0    K1    K2    K3    K3/4    K4
        -1.0    155    100    85    155    75    63    45    36    27
        -0.9    161    104    89    161    78    66    47    38    28
        -0.8    168    109    92    168    82    68    49    39    29
        -0.7    176    114    97    176    85    72    51    41    31
        -0.6    185    119    101    185    89    75    54    43    32
        -0.5    194    125    106    194    94    79    56    45    34
        -0.4    204    132    112    204    99    83    59    47    36
        -0.3    215    139    118    215    104    87    62    50    38
        -0.2    228    147    125    228    110    93    66    53    40
        -0.1    242    156    133    242    117    98    70    56    42
        0.0    258    167    142    258    125    105    75    60    45
        0.1    277    182    155    277    136    115    82    65    49
        0.2    298    200    170    298    150    126    90    72    54
        0.3    323    222    189    323    167    140    100    80    60
        0.4    352    250    212    352    187    157    112    90    67
        0.5    387    286    243    387    214    180    129    103    77
        0.6    431    333    283    431    250    210    150    120    90
        0.7    484    400    340    484    300    252    180    144    108
        0.8    554    500    425    554    375    315    225    180    135
        0.9    999    999    850    999    750    630    450    360    270
        1.0    999    999    999    999    999    999    999    999    999"""
        # self.Table = [[str(y) for y in x[8:].split("\t")] for x in self.String[1:].split("\n")]
        self.Table = [x[8:].split("\t") for x in self.String[1:].split("\n")]
        self.CraneGroup = dict(zip(["B0","B1","B2","B3","B4","B5","B6","B7"],range(2,21*7+2,21+2)))
        self.KFactor = dict(zip([str(x/10.0) for x in range(-10,11)],range(21)))
        self.NotchCase = dict(zip(["K","W0","W1","W2","K0","K1","K2","K3","K3/4","K4"],range(10)))
    def Get(self,C,K,N):
        return float(self.Table[self.CraneGroup[C]+self.KFactor[K]][self.NotchCase[N]])

FatigueStresses = cFatigueStresses()

def unittest():
    # getNormalDirection
    # Tensor
    # CalculateFatigue
    # FatigueStresses.Get
    pass

def destroyDB(result,step):
    ExtAPI.Log.WriteMessage("{0} Step: {3}, End Master Id: {1}, Item: {2}".format(str(datetime.now()), str(result.Properties.GetByName("Identifier").Value), str(result.Caption), str(step)))

class SlaveProperties:
    def __init__(self,result):
        Properties      = result.Properties
        self.Identifier      = Properties.GetByName("Identifier").Value
        self.Type            = Properties.GetByName("FatigueCalculation").Value
        self.Item            = Properties.GetByName("Item").Value
        self.WeldType        = Properties.GetByName("WeldType").Value
        self.NotchCase       = Properties.GetByName("NotchCase").Value
        self.StaticLoadCases = eval(Properties.GetByName("lc").Value)
        self.CraneGroup      = Properties.GetByName("CraneGroup").Value
        self.FatigueLoadCases = eval(Properties.GetByName("lcFatigue").Value)
    
    
def startSlave(result,step):
    # check if this is the last step, doe het alleen op de last step?
    global sp
    if result.Analysis.ResultsData.ResultSetCount != step:
        sp = None
        return
    ExtAPI.Log.WriteMessage("{0} Step: {3}, Start Slave Id: {1}, Item: {2}".format(str(datetime.now()), str(result.Properties.GetByName("Identifier").Value), str(result.Caption), str(step)))
    sp = SlaveProperties(result)
    
def endSlave(result,step):
    ExtAPI.Log.WriteMessage("{0} Step: {3}, End Slave Id: {1}, Item: {2}".format(str(datetime.now()), str(result.Properties.GetByName("Identifier").Value), str(result.Caption), str(step)))
    global sp
    del sp
    

# import doctest
# doctest.testmod()    
    