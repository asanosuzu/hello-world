# -*- coding: mbcs -*-
# -*- coding: utf-8 -*-
from numpy import *
import numpy as nmp



import os
import math

from math import *
from abaqus import *
from abaqusConstants import *
from caeModules import *
import numpy as np
from scipy import integrate
import sketch


def sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR):
    thick_list=[]
    thick_list1=[]
    thick_list2=[]
    r0=min(x_coordinate) 
    R=max(x_coordinate)
    alpha=arcsin(r0/R)*180.0/pi
    rb=r0+b_fiber
    r2b=r0+2*b_fiber    
    tR=2*nR*tP
    mR=2*pi*R*np.cos(alpha*pi/180.0)/b_fiber #赤道圆处缠绕纱片数
    #m0=2*pi*r0*np.cos(alpha*pi/180.0)/b_fiber #极孔处缠绕纱片数
    n0=1
    m0=mR*nR/n0

    int1,err1=integrate.quad(lambda r:2*pi*r*mR*nR/pi*np.arccos(r0/r)*tP,r0,rb)
    int2,err1=integrate.quad(lambda r:2*pi*r*mR*nR/pi*(np.arccos(r0/r2b)-np.arccos((r0+b_fiber)/r2b))*tP,rb,r2b)
    Vcon=int1+int2

    b1=tR*pi*R*np.cos(alpha*pi/180.0)/(m0*b_fiber)
    b2=mR*nR/pi*(np.arccos(r0/r2b)-np.arccos((r0+b_fiber)/r2b))*tP
    b3=mR*nR/pi*(r0/(r2b*(r2b**2-r0**2)**0.5)-rb/(r2b*(r2b**2-rb**2)**0.5))*tP
    b4=Vcon

    A=np.array([[1,r0,r0**2,r0**3],[1,r2b,r2b**2,r2b**3],[0,1,2*r2b,3*r2b**2],[pi*(r2b**2-r0**2),2*pi/3*(r2b**3-r0**3),pi/2*(r2b**4-r0**4),2*pi/5*(r2b**5-r0**5)]],dtype=float)
    B=np.array([[b1],[b2],[b3],[b4]],dtype=float)

    inv_A=np.linalg.inv(A)
    M=np.dot(inv_A,B)
    
    for i in x_coordinate:
        if i>=r0 and i<=r2b:
            partsss=i
            ti1=M[0]+M[1]*i+M[2]*i**2+M[3]*i**3
            thick_list1.append(float(ti1))

        elif i>=r2b and i<=R:
            ti2=mR*nR/pi*(np.arccos(r0/i)-np.arccos((r0+b_fiber)/i))*tP
            thick_list2.append(ti2)  
    thick_list=thick_list1+thick_list2
    
    new_thick_list=[i*(2*tP*nR)/thick_list[-1] for i in thick_list]
    #new_thick_list[0] = 0.01

    return new_thick_list,r0,R,alpha


def strengtheningThick(x_coordinate,b_fiber,tP,nR,Rs):
    thick_list=[]
    thick_list1=[]
    thick_list2=[]
    r0=min(x_coordinate) 
    #R=max(x_coordinate)
    alpha=arcsin(r0/Rs)*180.0/pi
    rb=r0+b_fiber
    r2b=r0+2*b_fiber    
    tR=2*nR*tP
    mR=2*pi*Rs*np.cos(alpha*pi/180.0)/b_fiber #赤道圆处缠绕纱片数
    #m0=2*pi*r0*np.cos(alpha*pi/180.0)/b_fiber #极孔处缠绕纱片数
    n0=1
    m0=mR*nR/n0

    int1,err1=integrate.quad(lambda r:2*pi*r*mR*nR/pi*np.arccos(r0/r)*tP,r0,rb)
    int2,err1=integrate.quad(lambda r:2*pi*r*mR*nR/pi*(np.arccos(r0/r2b)-np.arccos((r0+b_fiber)/r2b))*tP,rb,r2b)
    Vcon=int1+int2

    b1=tR*pi*Rs*np.cos(alpha*pi/180.0)/(m0*b_fiber)
    b2=mR*nR/pi*(np.arccos(r0/r2b)-np.arccos((r0+b_fiber)/r2b))*tP
    b3=mR*nR/pi*(r0/(r2b*(r2b**2-r0**2)**0.5)-rb/(r2b*(r2b**2-rb**2)**0.5))*tP
    b4=Vcon

    A=np.array([[1,r0,r0**2,r0**3],[1,r2b,r2b**2,r2b**3],[0,1,2*r2b,3*r2b**2],[pi*(r2b**2-r0**2),2*pi/3*(r2b**3-r0**3),pi/2*(r2b**4-r0**4),2*pi/5*(r2b**5-r0**5)]],dtype=float)
    B=np.array([[b1],[b2],[b3],[b4]],dtype=float)

    inv_A=np.linalg.inv(A)
    M=np.dot(inv_A,B)
    
    for i in x_coordinate:
        if i>=r0 and i<=r2b:
            partsss=i
            ti1=M[0]+M[1]*i+M[2]*i**2+M[3]*i**3
            thick_list1.append(float(ti1))

        else:#if i>=r2b and i<=Rs:
            ti2=mR*nR/pi*(np.arccos(r0/i)-np.arccos((r0+b_fiber)/i))*tP
            thick_list2.append(ti2)  
    thick_list=thick_list1+thick_list2
    
    new_thick_list=[i*(2*tP*nR)/thick_list[-1] for i in thick_list]
 

    return new_thick_list,r0,alpha




def mian(L_style,slopeDegree,segmentationNum,rotationAngle,fibersWidth,
        pieceNumber,approximateSize,typeReinforce,
        strengtheningAngle,modelName,CompositeName,strThick,assemblyName,
        L_setName,R_setName):

    L_tables=[]
    R_tables=[]
    L_helixThickness,R_helixThickness=[],[]
    L_hoopThickness,R_hoopThickness=[],[]
    L_strengthThickness,R_strengthThickness=[],[]
    L_reamingThickness,R_reamingThickness=[],[]
    L_reamingAngle=[]
    R_reamingAngle=[]
    L_strengthLocation=[]
    R_strengthLocation=[]
    
    for x in L_style:
        if x[5] == "up":
            if x[0] == 'helix':
                L_tables.append('helix')
                L_helixThickness.append(x[4])
                #L_reamingAngle.append(x[1])
            elif x[0] == 'helixReaming':
                L_tables.append('helixReaming')
                L_reamingAngle.append(x[1])
                L_reamingThickness.append(x[4])
            elif x[0] == 'strengthening':
                L_tables.append('strengthening')
                L_strengthLocation.append((x[2],x[3]))
                L_strengthThickness.append(x[4])
            elif x[0] == 'hoop':
                L_tables.append('hoop')
                L_hoopThickness.append(x[4])
        else:
            if x[0] == 'helix':
                R_tables.append('helix')
                R_helixThickness.append(x[4])
                #L_reamingAngle.append(x[1])
            elif x[0] == 'helixReaming':
                R_tables.append('helixReaming')
                R_reamingAngle.append(x[1])
                R_reamingThickness.append(x[4])
            elif x[0] == 'strengthening':
                R_tables.append('strengthening')
                R_strengthLocation.append((x[2],x[3]))
                R_strengthThickness.append(x[4])
            elif x[0] == 'hoop':
                R_tables.append('hoop')
                R_hoopThickness.append(x[4])
            
    
    
        

    xiedegree=slopeDegree #环向层斜度<
    todegree=pi/180.0

    huanxiang_alpha=90.0 #hoopAlpha #环向角度
    num=segmentationNum#（封头切片数）不加.0
    rotangle=rotationAngle#旋转角度
     
    b_fiber=fibersWidth #纱线宽度
    # tP=singleFiberThick #单层纱带的厚度
    # siglethick=hoopSingleThick#环向层单层厚度
    nR=1 #赤道处螺旋缠绕层数


    Error_control=1e-6
    piece_number=pieceNumber
    approximate_size=approximateSize 

    if typeReinforce =="windReinforce":
        type_reinforce="chanraobuqiang"##"pufangbuqiang" #
    else:
        type_reinforce="pufangbuqiang"
        
    buqiangangle=strengtheningAngle
    b_thick=strThick
    modelName = modelName
    #linerName = LinerModel
    compositeName = CompositeName
    m = mdb.models[modelName]

    #创建复合材料工程常数
    m.HomogeneousSolidSection(name='Section-composite', material=compositeName, thickness=None)
    #m.HomogeneousSolidSection(name='Section-liner', material=linerName, thickness=None)

	


    curspline=[]
    curve_x=[]
    curve_y=[]
    toppoint=[]
    botpoint=[]
    all_angle=[]
    all_count=[]
    all_name=[]
    part_name=[]
    all_instances=[]
    tangentpoint={}
    a=m.rootAssembly
    i=a.instances[assemblyName]
    edgeset = a.sets[L_setName]
    #outlineset=a.sets['Set-2']
    ver=i.vertices
    L_liner = m.ConstrainedSketch(name='L_liner', sheetSize=200.0)
    g, v, d, c = L_liner.geometry, L_liner.vertices, L_liner.dimensions, L_liner.constraints
    cline=L_liner.ConstructionLine((0,10),(0,-10))
    L_liner.assignCenterline(line=cline)
    ###读取内衬数据

    liner_line=[]
    all_spline=[]
    edpoint=[]

    for edge in edgeset.edges:
        try:
            lineLength=edge.getSize(printResults=False)
            partnum=num
            for i in range(1,partnum):
                tonum=1.0/partnum*i
                curpoint = edge.getCurvature(parameter=tonum)
                evaPoint=curpoint['evaluationPoint']
                tabgent=curpoint['tangent']
                tangentpoint[(evaPoint[0],evaPoint[1])]=tabgent
                curve_x.append(evaPoint[0])
                curve_y.append(evaPoint[1])
            curveLineVertices=edge.getVertices()
            topindex,botindex=curveLineVertices[0],curveLineVertices[1]
            toppoint=ver[topindex].pointOn[0]
            botpoint=ver[botindex].pointOn[0]
            spline=zip(curve_x,curve_y) 
            spline.insert(0,(toppoint[0],toppoint[1]))
            spline.insert(-1,(botpoint[0],botpoint[1]))
            tangentpoint[(toppoint[0],toppoint[1])]=(0,1,1)
            tangentpoint[(botpoint[0],botpoint[1])]=(0,1,2)
            spline.sort(key=(lambda x: x[1]))
                  
            
            L_liner.Spline(points=spline) 
            curve_x=[]
            curve_y=[]
            spline=[]
        except:
            straightLineVertices=edge.getVertices()
            topindex,botindex=straightLineVertices[0],straightLineVertices[1]
            toppoint=ver[topindex].pointOn[0]
            botpoint=ver[botindex].pointOn[0]
            linecoord=[(toppoint[0]+botpoint[0])/2.0 ,(toppoint[1] + botpoint[1])/2.0]
            line1=L_liner.Line(point1=(toppoint[0], toppoint[1]), point2=(botpoint[0], botpoint[1])) 


    #画缠绕层
    x_composite=[]
    y_composite=[]
    x_coordinate=[]
    y_coordinate=[]


    incoordinate=tangentpoint.keys()
    incoordinate.sort(key=(lambda x: x[1]))
    incoordinate.reverse()
    for i in range(len(incoordinate)):
        x_coordinate.append(incoordinate[i][0])
        y_coordinate.append(incoordinate[i][1])
        
    compositei=zip(x_coordinate,y_coordinate)


    L = min(y_coordinate) #筒身长度一半
    R0 = max(x_coordinate)

    #r#扩孔半径
    siglethick = 0.2
    cankaolunkuo=zip(x_coordinate,y_coordinate)
    solvek=zip(x_coordinate,y_coordinate)
    y2=L-xiedegree
    #y2=L-siglethick/tan(xiedegree*todegree)
    leftAcc=y2
    #thicklist,r0,R0,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)
    #luoxuanhoudu=thicklist[-1]
    
    L_hxcount=0   
    L_lxcount=0 
    L_reaming=0 
    L_buqiang=0
    ftBotFaceArray=[]
    ftTopFaceArray=[]
    tsBotFaceArray=[]
    tsTopFaceArray=[]


    L_helixThicknessCount = 0
    L_reamingThicknessCount = 0
    L_strengthThicknessCount = 0
    L_hoopThicknessCount = 0
    
    L_thickness=R0

    zonghoudu=[]
    zongjiaodu=[]


    for windtype in L_tables:
        if windtype == 'helix':
            fengtoui=cankaolunkuo#内轮廓等于参考轮廓
            #更新实际坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1])
            tP = L_helixThickness[L_helixThicknessCount]
            L_helixThicknessCount+=1
            thick_list,r0,R,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)
            
            zonghoudu.append(np.array(tP))
            

            
            s1 = m.ConstrainedSketch(name='L_helixLayer'+str(L_lxcount-L_buqiang), sheetSize=200.0)

            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]
            for i in range(len(fengtoui)-1):
                k1=(solvek[i+1][1]-solvek[i][1])/(solvek[i+1][0]-solvek[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            xvalue.append(fengtoui[-1][0]+thick_list[-1])
            yvalue.append(fengtoui[-1][1])
            fengtouo=zip(xvalue,yvalue)
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
                

            
            #筒身
            newRi=L_thickness
            L_thickness+=tP*2
            newRo=newRi+tP*2
            line1=s1.Line(point1=(newRi,y2), point2=(newRi,0))
            line2=s1.Line(point1=(newRo,y2), point2=(newRo,0))
            #连线
            line1=s1.Line(point1=(newRi,y2), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            line2=s1.Line(point1=(newRo,y2), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line4=s1.Line(point1=(newRi,0), point2=(newRo,0))
            L_lxcount=L_lxcount+1
            #更新参考轮廓
            cankaolunkuo=fengtouo
            p1 = m.Part(name='L_helixLayer'+str(L_lxcount-L_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
                   
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='L_helixLayer'+str(L_lxcount-L_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('L_helixLayer'+str(L_lxcount-L_buqiang)+'-1') 

            #设置材料赋予参考面
            f=p1.faces
            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='L_surfHeadTop'+str(L_lxcount)+'-'+str(i+1))
            side1Faces22 = f.findAt(((newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces22, name='L_surfCylinderTop'+str(L_lxcount))        
            #封头画切分线
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickface=f1.findAt(coordinates=(0.5*(newRi+newRo), y2/2.0, 0.0))
            pickedge=e1.findAt(coordinates=(newRo, y2/2.0, 0.0))
            t = p1.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.5*(newRi+newRo), 
                y2/2.0, 0.0)), sketchUpEdge=e1.findAt(coordinates=(newRo, 
                y2/2.0, 0.0)), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=650.45, gridSpacing=16.26, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            for k in range(1,len(fengtoui)):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=pickedge, faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']   
            #封头处切分
            p1 = m.parts['L_helixLayer'+str(L_lxcount-L_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度

            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            

            
            for k in face_coordinate:            
                count=count+1
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                winding_angle=arcsin(r0/k[0])*180.0/pi
                list_angle.append(winding_angle)
                p1.Set(name='L_set-'+str(L_lxcount)+'-'+str(count),cells=pickcells)
            
            all_angle.append(list_angle)
            all_count.append(count)  
          
 
            SDcoordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))

        elif windtype == 'helixReaming':
            kfengtoui=cankaolunkuo#内轮廓等于参考轮廓
            reamingangle_1=L_reamingAngle[L_reaming]
            L_reaming+=1
            temp_x = []
            for i in range(len(kfengtoui)):
                temp_x.append(kfengtoui[i][0])
            Rs = max(temp_x)
            r=Rs*sin(reamingangle_1*pi/180)

            indexnum = 0
            #判定扩孔的位置
            for i in range(len(kfengtoui)-1):
                a=kfengtoui[i][0]
                b=kfengtoui[i+1][0]
                if a<=r and b>r:
                    indexnum=i
            if indexnum == 0:
                msg = "Check the reaming Angle" 
                raise AbaqusException, msg
                
                
    # tP=singleFiberThick #单层纱带的厚度
    # siglethick=hoopSingleThick#环向层单层厚度                
                
            #扩孔的实际内轮廓
            fengtoui=kfengtoui[indexnum:]
            #右边尖孔
            left_point = kfengtoui[indexnum-1]
            solvek1=solvek[indexnum:]
            #更新实际位置坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1]) 
            
            tP = L_reamingThickness[L_reamingThicknessCount]
            L_reamingThicknessCount+=1
            thick_list,r0,R,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)
            s1 = m.ConstrainedSketch(name='L_helixLayer'+str(L_lxcount-L_buqiang), sheetSize=200.0)
            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)

     
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]
            for i in range(len(fengtoui)-1):
                k1=(solvek1[i+1][1]-solvek1[i][1])/(solvek1[i+1][0]-solvek1[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            xvalue.append(fengtoui[-1][0]+thick_list[-1])
            yvalue.append(fengtoui[-1][1])
            fengtouo=zip(xvalue,yvalue)
            
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
            
            #s1.Spline(points=fengtoui)    
            #s1.Spline(points=fengtouo)      
            
            #筒身
            newRi=L_thickness
            L_thickness+=tP*2
            newRo=newRi+tP*2
            line1=s1.Line(point1=(newRi,y2), point2=(newRi,0))
            line2=s1.Line(point1=(newRo,y2), point2=(newRo,0))
            #封头连线
            line1=s1.Line(point1=(newRi,y2), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            line2=s1.Line(point1=(newRo,y2), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            #line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line3_1 = s1.Line(point1=(fengtouo[0][0],fengtouo[0][1]),point2=(left_point[0],left_point[1]))
            line3_2 = s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]),point2=(left_point[0],left_point[1]))
            line4=s1.Line(point1=(newRi,0), point2=(newRo,0))
            L_lxcount=L_lxcount+1
            #更新参考轮廓
            cankaolunkuo=kfengtoui[:indexnum]+fengtouo
            p1 = m.Part(name='L_helixLayer'+str(L_lxcount-L_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='L_helixLayer'+str(L_lxcount-L_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('L_helixLayer'+str(L_lxcount-L_buqiang)+'-1')

            #设置材料赋予参考面
            f=p1.faces

            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='L_surfHeadTop'+str(L_lxcount)+'-'+str(i+1))        
            

            side1Faces22 = f.findAt(((newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces22, name='L_surfCylinderTop'+str(L_lxcount))        
            #封头画切分线
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickface=f1.findAt(coordinates=(0.5*(newRi+newRo), y2/2.0, 0.0))
            pickedge=e1.findAt(coordinates=(newRo, y2/2.0, 0.0))
            t = p1.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.5*(newRi+newRo), 
                y2/2.0, 0.0)), sketchUpEdge=e1.findAt(coordinates=(newRo, 
                y2/2.0, 0.0)), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=650.45, gridSpacing=16.26, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            for k in range(1,len(fengtoui)):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=pickedge, faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']
            #封头处切分
            p1 = m.parts['L_helixLayer'+str(L_lxcount-L_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度

            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            for k in face_coordinate:
                count=count+1
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                #print r,r0,(k[0],k[1], k[2])
                winding_angle=arcsin(r0/k[0])*180.0/pi
                list_angle.append(winding_angle)
                p1.Set(name='L_set-'+str(L_lxcount)+'-'+str(count),cells=pickcells)
            all_angle.append(list_angle)    
            #赋予角度
            all_count.append(count)
            
            
            #切分斜角多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[0], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[0][0], 
                fengtoui[0][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[0][0], fengtouo[0][1], 
                0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[0], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[0][0]+fengtouo[0][0]), 0.5*(fengtoui[0][1]+fengtouo[0][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), 
                L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=pickedCells, edges=pickedEdges)            
           
            
            SDcoordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))
            # fcoordinate=side1Faces22.pointsOn[0][0]
            # c = p1.cells
            # f = p1.faces
            # p1.assignStackDirection(referenceRegion=f.findAt(coordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))), cells=c)
        


        elif windtype == 'strengthening':
            kfengtoui=cankaolunkuo#内轮廓等于参考轮廓
            temp_x = []
            for i in range(len(kfengtoui)):
                temp_x.append(kfengtoui[i][0])
            Rs = max(temp_x)
            x1=L_strengthLocation[L_buqiang][0]
            x2=L_strengthLocation[L_buqiang][1]
            
            indexnum_i,indexnum_j =0,0
            
            if x2 <= L:
                indexnum_j=len(kfengtoui)
                # right_point = kfengtoui[-1]
                
            else:
                for j in range(len(kfengtoui)-1):
                    c=kfengtoui[j][1]
                    d=kfengtoui[j+1][1]
                    if c>=x2 and d<x2:
                        indexnum_j=j
                    
                right_point = kfengtoui[indexnum_j+1]
                
            #判定补强左位置
            
            for i in range(len(kfengtoui)-1):
                a=kfengtoui[i][1]
                b=kfengtoui[i+1][1]
                if a>=x1 and b<x1:
                    indexnum_i=i
                    
             #判定补强右位置

            if indexnum_i == 0 or indexnum_j == 0:
                msg = "Check the strengthening location" 
                raise AbaqusException, msg
                        

            
                
            #扩孔的实际内轮廓
            fengtoui=kfengtoui[indexnum_i:indexnum_j+1]
            #左边和右边尖孔
            left_point = kfengtoui[indexnum_i-1]
            #right_point = kfengtoui[indexnum_j-1]
            r0=left_point[0]
            solvek1=solvek[indexnum_i:indexnum_j+1]

            #print len(solvek1)
            #更新实际位置坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1])     
            if type_reinforce == "chanraobuqiang":
                #print len(x_coordinate)
                tP = L_strengthThickness[L_strengthThicknessCount]
                L_strengthThicknessCount+=1
                thick_list,r0,alpha= strengtheningThick(x_coordinate,b_fiber,tP,nR,Rs)   
                
            else:
                thick_list = [b_thick]*len(fengtoui)
            
            L_buqiang=L_buqiang+1
            s1 = m.ConstrainedSketch(name='L_strengtheningLayer'+str(L_buqiang), sheetSize=200.0)
            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)

     
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]
            for i in range(len(fengtoui)-1):
                k1=(solvek1[i+1][1]-solvek1[i][1])/(solvek1[i+1][0]-solvek1[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            # xvalue.append(fengtoui[-1][0]+thick_list[-1]*cos(theta1))
            # yvalue.append(fengtoui[-1][0]+thick_list[-1]*sin(theta1))
            fengtouo=zip(xvalue,yvalue)
            
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
            
            #s1.Spline(points=fengtoui)    
            #s1.Spline(points=fengtouo)      
            
            #连线
            # line1=s1.Line(point1=(right_point[0],right_point[1]), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            # line2=s1.Line(point1=(right_point[0],right_point[1]), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            line2=s1.Line(point1=(fengtoui[-1][0],fengtoui[-1][1]), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            #line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line3_1 = s1.Line(point1=(fengtouo[0][0],fengtouo[0][1]),point2=(left_point[0],left_point[1]))
            line3_2 = s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]),point2=(left_point[0],left_point[1]))

            
            L_lxcount=L_lxcount+1
            #更新参考轮廓
            if x2 <= L:
                cankaolunkuo1=kfengtoui[:indexnum_i]+fengtouo
                cankaolunkuo1.append(fengtoui[-1])
                cankaolunkuo=cankaolunkuo1+kfengtoui[indexnum_j:]
            else:
                cankaolunkuo=kfengtoui[:indexnum_i]+fengtouo+kfengtoui[indexnum_j:]
            
            p1 = m.Part(name='L_strengtheningLayer'+str(L_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='L_strengtheningLayer'+str(L_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('L_strengtheningLayer'+str(L_buqiang)+'-1')

            #存储侧边参考面
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.pointOn
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            #设置材料赋予参考面
            f=p1.faces

            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='L_surfStrengtheningTop'+str(L_lxcount)+'-'+str(i+1))        
            

            
             
            p1 = m.parts['L_strengtheningLayer'+str(L_buqiang)]
            plant1=p1.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
            d1=p1.datums
            
            #封头画切分线
            t = p1.MakeSketchTransform(sketchPlane=f.findAt(coordinates=centercoordinate[0]
            ), sketchUpEdge=d1[plant1.id], sketchPlaneSide=SIDE1, origin=(0.0, 
                0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=1258.72, gridSpacing=31.46, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
         
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            
            pickface=f.findAt(coordinates=centercoordinate[0])
            
            
            for k in range(1,len(fengtoui)-2):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=d1[plant1.id], faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']            
            #封头处切分
            p1 = m.parts['L_strengtheningLayer'+str(L_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)-2):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度
            
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            c=p1.cells
            if type_reinforce == "chanraobuqiang":
                for k in face_coordinate:
                    count=count+1
                    pickcells=c.findAt(((k[0], k[1], k[2]),))
                    #print "r0=%s,k[0]=%s\n"%(r0,k[0])
                    winding_angle=arcsin(r0/k[0])*180.0/pi
                    list_angle.append(winding_angle)
                    p1.Set(name='L_set-'+str(L_lxcount)+'-'+str(count),cells=pickcells)
                all_angle.append(list_angle)   
            else:
                for k in face_coordinate:
                    count=count+1
                    pickcells=c.findAt(((k[0], k[1], k[2]),))
                    winding_angle=buqiangangle

                    list_angle.append(winding_angle)
                    p1.Set(name='L_set-'+str(L_lxcount)+'-'+str(count),cells=pickcells)
                all_angle.append(list_angle)   
            #赋予角度
            all_count.append(count)


            #切分斜角前面多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[0], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[0][0], 
                fengtoui[0][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[0][0], fengtouo[0][1], 0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[0], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[0][0]+fengtouo[0][0]), 0.5*(fengtoui[0][1]+fengtouo[0][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], 
                -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges) 


    #        切分斜角后面多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[-1], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[-2][0], 
                fengtoui[-2][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[-1][0], fengtouo[-1][1], 0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[-1], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[-2][0]+fengtouo[-1][0]), 0.5*(fengtoui[-2][1]+fengtouo[-1][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], 
                -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges) 


            
 

        

        else:           
            s2 = m.ConstrainedSketch(name='L_hoopLayer'+str(L_hxcount), sheetSize=200.0)
            cline=s2.ConstructionLine((0,10),(0,-10))
            #筒身

            tP = L_hoopThickness[L_hoopThicknessCount]
            
            tnewRi=L_thickness
            L_thickness+=tP
            tnewRo=tnewRi+tP

            line1=s2.Line(point1=(tnewRi,y2), point2=(tnewRi,0))
            line2=s2.Line(point1=(tnewRo,y2), point2=(tnewRo,0))        
            #连线
            line1=s2.Line(point1=(cankaolunkuo[-1][0],cankaolunkuo[-1][1]), point2=(tnewRi,y2))
            line2=s2.Line(point1=(cankaolunkuo[-1][0],cankaolunkuo[-1][1]), point2=(tnewRo,y2))
            line3=s2.Line(point1=(tnewRi,0), point2=(tnewRo,0))
            L_hxcount=L_hxcount+1
            p1 = m.Part(name='L_hoopLayer'+str(L_hxcount), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s2, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)        
            #设置材料赋予参考面
            f=p1.faces
            side1Faces1 = f.findAt(((tnewRi*cos(0.5*rotangle*todegree),y2/2.0,-tnewRi*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces1, name='L_surfHoopBottom')       
            side1Faces2 = f.findAt(((tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces2, name='L_surfHoopTop'+str(L_hxcount))
            #设置装配绑定参考面 
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            huanxiangceng=a1.Instance(name='L_hoopLayer'+str(L_hxcount)+'-1', part=p1, dependent=ON)  
            all_instances.append(huanxiangceng)
            all_name.append('L_hoopLayer'+str(L_hxcount)+'-1')
            c=p1.cells
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.pointOn
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0]) 
            for k in face_coordinate:
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                p1.Set(name='L_set-'+str(L_hxcount),cells=pickcells)

            SDcoordinates=(tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree))
            # c = p1.cells
            # f = p1.faces
            # p1.assignStackDirection(referenceRegion=f.findAt(coordinates=(tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree))), cells=c)


    ###############################
    ###########################################
    #########################################################
    #################################################################

    #all_angle
    curspline=[]
    curve_x=[]
    curve_y=[]
    toppoint=[]
    botpoint=[]
    r_all_angle=[]
    r_all_count=[]
    all_name=[]
    #all_instances=[]
    tangentpoint={}
    a=m.rootAssembly
    i=a.instances[assemblyName]
    edgeset = a.sets[R_setName]
    #outlineset=a.sets['Set-2']
    ver=i.vertices
    s_liner = m.ConstrainedSketch(name='R_liner1', sheetSize=200.0)
    g, v, d, c = s_liner.geometry, s_liner.vertices, s_liner.dimensions, s_liner.constraints
    cline=s_liner.ConstructionLine((0,10),(0,-10))
    s_liner.assignCenterline(line=cline)
    ###读取内衬数据

    liner_line=[]
    all_spline=[]
    edpoint=[]




    for edge in edgeset.edges:
        try:
            lineLength=edge.getSize(printResults=False)
            if lineLength<=20.0:
                partnum=20
            else:
                partnum=num
            for i in range(1,partnum):
                tonum=1.0/partnum*i
                curpoint = edge.getCurvature(parameter=tonum)
                evaPoint=curpoint['evaluationPoint']
                tabgent=curpoint['tangent']
                tangentpoint[(evaPoint[0],evaPoint[1])]=tabgent
                curve_x.append(evaPoint[0])
                curve_y.append(evaPoint[1])
            curveLineVertices=edge.getVertices()
            topindex,botindex=curveLineVertices[0],curveLineVertices[1]
            toppoint=ver[topindex].pointOn[0]
            botpoint=ver[botindex].pointOn[0]
            spline=zip(curve_x,curve_y) 
            spline.insert(0,(toppoint[0],toppoint[1]))
            spline.insert(-1,(botpoint[0],botpoint[1]))
            tangentpoint[(toppoint[0],toppoint[1])]=(0,1,1)
            tangentpoint[(botpoint[0],botpoint[1])]=(0,1,2)
            spline.sort(key=(lambda x: x[1]))
                  
            
            #s_liner.Spline(points=spline) 
            curve_x=[]
            curve_y=[]
            spline=[]
        except:
            straightLineVertices=edge.getVertices()
            topindex,botindex=straightLineVertices[0],straightLineVertices[1]
            toppoint=ver[topindex].pointOn[0]
            botpoint=ver[botindex].pointOn[0]
            linecoord=[(toppoint[0]+botpoint[0])/2.0 ,(toppoint[1] + botpoint[1])/2.0]
            #line1=s_liner.Line(point1=(toppoint[0], toppoint[1]), point2=(botpoint[0], botpoint[1])) 








    #画缠绕层

    x_composite=[]
    y_composite=[]
    x_coordinate=[]
    y_coordinate=[]


     

    incoordinate=tangentpoint.keys()
    incoordinate.sort(key=(lambda x: x[1]))
    #incoordinate.reverse()
    for i in range(len(incoordinate)):
        x_coordinate.append(incoordinate[i][0])
        y_coordinate.append(incoordinate[i][1])
        
    compositei=zip(x_coordinate,y_coordinate)



    L=max(y_coordinate) #筒身长度一半



    #r#扩孔半径


    cankaolunkuo=zip(x_coordinate,y_coordinate)
    solvek=zip(x_coordinate,y_coordinate)
    y2=L+xiedegree
    #y2=L+siglethick/tan(xiedegree*todegree)
    rightAcc=y2
    thicklist,r0,R0,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)
    luoxuanhoudu=thicklist[-1]
      
    R_hxcount=0   
    R_lxcount=0 
    R_buqiang=0
    R_reaming=0
    ftBotFaceArray=[]
    ftTopFaceArray=[]
    tsBotFaceArray=[]
    tsTopFaceArray=[]


    R_helixThicknessCount = 0
    R_reamingThicknessCount = 0
    R_strengthThicknessCount = 0
    R_hoopThicknessCount = 0
    
    R_thickness=R0



    zonghoudu=[]
    zongjiaodu=[]
    #b_thick=0.3

    for windtype in R_tables:
        if windtype == 'helix':
            fengtoui=cankaolunkuo#内轮廓等于参考轮廓
            #更新实际坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1])
            tP = R_helixThickness[R_helixThicknessCount]
            R_helixThicknessCount+=1
            thick_list,r0,R,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)
            zonghoudu.append(np.array(tP))
            s1 = m.ConstrainedSketch(name='R_helixLayer'+str(R_lxcount-R_buqiang), sheetSize=200.0)

            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]
            for i in range(len(fengtoui)-1):
                k1=(solvek[i+1][1]-solvek[i][1])/(solvek[i+1][0]-solvek[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            xvalue.append(fengtoui[-1][0]+thick_list[-1])
            yvalue.append(fengtoui[-1][1])
            fengtouo=zip(xvalue,yvalue)
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
                
            #s1.Spline(points=fengtoui)    
            #s1.Spline(points=fengtouo)      
            
            #筒身
            newRi=R_thickness
            R_thickness+=tP*2
            newRo=newRi+tP*2
            line1=s1.Line(point1=(newRi,y2), point2=(newRi,0))
            line2=s1.Line(point1=(newRo,y2), point2=(newRo,0))
            #连线
            line1=s1.Line(point1=(newRi,y2), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            line2=s1.Line(point1=(newRo,y2), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line4=s1.Line(point1=(newRi,0), point2=(newRo,0))
            R_lxcount=R_lxcount+1
            #更新参考轮廓
            cankaolunkuo=fengtouo
            p1 = m.Part(name='R_helixLayer'+str(R_lxcount-R_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
                   
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='R_helixLayer'+str(R_lxcount-R_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('R_helixLayer'+str(R_lxcount-R_buqiang)+'-1') 

            #设置材料赋予参考面
            f=p1.faces
            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='R_surfHeadTop'+str(R_lxcount)+'-'+str(i+1))
            side1Faces22 = f.findAt(((newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces22, name='R_surfCylinderTop'+str(R_lxcount))        
            #封头画切分线
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickface=f1.findAt(coordinates=(0.5*(newRi+newRo), y2/2.0, 0.0))
            pickedge=e1.findAt(coordinates=(newRo, y2/2.0, 0.0))
            t = p1.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.5*(newRi+newRo), 
                y2/2.0, 0.0)), sketchUpEdge=e1.findAt(coordinates=(newRo, 
                y2/2.0, 0.0)), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=650.45, gridSpacing=16.26, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            for k in range(1,len(fengtoui)):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=pickedge, faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']            
            #封头处切分
            p1 = m.parts['R_helixLayer'+str(R_lxcount-R_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度

            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            for k in face_coordinate:            
                count=count+1
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                winding_angle=arcsin(r0/k[0])*180.0/pi
                list_angle.append(winding_angle)
                p1.Set(name='R_set-'+str(R_lxcount)+'-'+str(count),cells=pickcells)
            
            r_all_angle.append(list_angle)
            r_all_count.append(count)    

            SDcoordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))
            # fcoordinate=side1Faces22.pointsOn[0][0]
            # c = p1.cells
            # f = p1.faces
            # p1.assignStackDirection(referenceRegion=f.findAt(coordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))), cells=c)     

        elif windtype == 'helixReaming':
            kfengtoui=cankaolunkuo#内轮廓等于参考轮廓
            reamingangle_1=R_reamingAngle[R_reaming]
            R_reaming+=1
            temp_x = []
            for i in range(len(kfengtoui)):
                temp_x.append(kfengtoui[i][0])
            Rs = max(temp_x)
            
            r=Rs*sin(reamingangle_1*pi/180)
            indexnum = 0
            #判定扩孔的位置
            for i in range(len(kfengtoui)-1):
                a=kfengtoui[i][0]
                b=kfengtoui[i+1][0]
                if a<=r and b>r:
                    indexnum=i
                    
            if indexnum == 0:
                msg = "Check the reaming Angle" 
                raise AbaqusException, msg                    
            #扩孔的实际内轮廓
            fengtoui=kfengtoui[indexnum:]
            #右边尖孔
            left_point = kfengtoui[indexnum-1]
            solvek1=solvek[indexnum:]
            #更新实际位置坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1])     
            
            tP = R_reamingThickness[R_reamingThicknessCount]
            R_reamingThicknessCount+=1
            thick_list,r0,R,alpha=sanciyangtiaohouduyuce(x_coordinate,b_fiber,tP,nR)

            s1 = m.ConstrainedSketch(name='R_helixLayer'+str(R_lxcount-R_buqiang), sheetSize=200.0)
            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)

     
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]
            for i in range(len(fengtoui)-1):
                k1=(solvek1[i+1][1]-solvek1[i][1])/(solvek1[i+1][0]-solvek1[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            xvalue.append(fengtoui[-1][0]+thick_list[-1])
            yvalue.append(fengtoui[-1][1])
            fengtouo=zip(xvalue,yvalue)
            
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
            
            #s1.Spline(points=fengtoui)    
            #s1.Spline(points=fengtouo)      
            
            #筒身
            newRi=R_thickness
            R_thickness+=tP*2
            newRo=newRi+tP*2
            line1=s1.Line(point1=(newRi,y2), point2=(newRi,0))
            line2=s1.Line(point1=(newRo,y2), point2=(newRo,0))
            #连线
            line1=s1.Line(point1=(newRi,y2), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            line2=s1.Line(point1=(newRo,y2), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            #line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line3_1 = s1.Line(point1=(fengtouo[0][0],fengtouo[0][1]),point2=(left_point[0],left_point[1]))
            line3_2 = s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]),point2=(left_point[0],left_point[1]))
            line4=s1.Line(point1=(newRi,0), point2=(newRo,0))
            R_lxcount=R_lxcount+1
            #更新参考轮廓
            cankaolunkuo=kfengtoui[:indexnum]+fengtouo
            p1 = m.Part(name='R_helixLayer'+str(R_lxcount-R_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='R_helixLayer'+str(R_lxcount-R_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('R_helixLayer'+str(R_lxcount-R_buqiang)+'-1')

            #设置材料赋予参考面
            f=p1.faces

            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='R_surfHeadTop'+str(R_lxcount)+'-'+str(i+1))        
            

            side1Faces22 = f.findAt(((newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces22, name='R_surfCylinderTop'+str(R_lxcount))        
            #封头画切分线
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickface=f1.findAt(coordinates=(0.5*(newRi+newRo), y2/2.0, 0.0))
            pickedge=e1.findAt(coordinates=(newRo, y2/2.0, 0.0))
            t = p1.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.5*(newRi+newRo), 
                y2/2.0, 0.0)), sketchUpEdge=e1.findAt(coordinates=(newRo, 
                y2/2.0, 0.0)), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=650.45, gridSpacing=16.26, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            for k in range(1,len(fengtoui)):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=pickedge, faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']            
            #封头处切分
            p1 = m.parts['R_helixLayer'+str(R_lxcount-R_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), 
                    L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度

            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            for k in face_coordinate:
                count=count+1
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                #print r,r0,(k[0],k[1], k[2])
                winding_angle=arcsin(r0/k[0])*180.0/pi
                list_angle.append(winding_angle)
                p1.Set(name='R_set-'+str(R_lxcount)+'-'+str(count),cells=pickcells)
            r_all_angle.append(list_angle)   
    #all_count        
            #赋予角度
            r_all_count.append(count)
            
            #切分斜角多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[0], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[0][0], 
                fengtoui[0][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[0][0], fengtouo[0][1], 
                0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[0], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[0][0]+fengtouo[0][0]), 0.5*(fengtoui[0][1]+fengtouo[0][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[-1][0]*cos(0.5*rotangle*todegree), 
                L, -1*fengtouo[-1][0]*sin(0.5*rotangle*todegree))), cells=pickedCells, edges=pickedEdges)         
             
            SDcoordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))
            # fcoordinate=side1Faces22.pointsOn[0][0]
            # c = p1.cells
            # f = p1.faces
            # p1.assignStackDirection(referenceRegion=f.findAt(coordinates=(newRo*cos(0.5*rotangle*todegree),y2/2.0,-newRo*sin(0.5*rotangle*todegree))), cells=c)
        


        elif windtype == 'strengthening':
            kfengtoui=cankaolunkuo#内轮廓等于参考轮廓
            x1=R_strengthLocation[R_buqiang][0]
            x2=R_strengthLocation[R_buqiang][1]
            
            temp_x = []
            for i in range(len(kfengtoui)):
                temp_x.append(kfengtoui[i][0])
            Rs = max(temp_x)            
            
            indexnum_j,indexnum_i =0,0
            if x2 >= L:
                indexnum_j=len(kfengtoui)
                # right_point = kfengtoui[-1]
                
            else:
                for j in range(len(kfengtoui)-1):
                    d=kfengtoui[j][1]
                    c=kfengtoui[j+1][1]
                    if c>=x2 and d<x2:
                        indexnum_j=j
                    
                #right_point = kfengtoui[indexnum_j+1]
                
            #判定补强左位置
            
            for i in range(len(kfengtoui)-1):
                b=kfengtoui[i][1]
                a=kfengtoui[i+1][1]
                if a>=x1 and b<x1:
                    indexnum_i=i
                    
            if indexnum_i == 0 or indexnum_j == 0:
                msg = "Check the strengthening location" 
                raise AbaqusException, msg
            
             #判定补强右位置
   
            #扩孔的实际内轮廓
            fengtoui=kfengtoui[indexnum_i:indexnum_j+1]
            #左边和右边尖孔
            left_point = kfengtoui[indexnum_i-1]
            #right_point = kfengtoui[indexnum_j-1]
            r0=left_point[0]
            solvek1=solvek[indexnum_i:indexnum_j+1]

            #print len(solvek1)
            #更新实际位置坐标
            x_coordinate=[]
            y_coordinate=[]
            for i in range(len(fengtoui)):
                x_coordinate.append(fengtoui[i][0])
                y_coordinate.append(fengtoui[i][1])     
            if type_reinforce == "chanraobuqiang":
                tP = R_strengthThickness[R_strengthThicknessCount]
                R_strengthThicknessCount+=1
                thick_list,r0,alpha= strengtheningThick(x_coordinate,b_fiber,tP,nR,Rs)   
                
            else:
                thick_list = [b_thick]*len(fengtoui)
            
            R_buqiang=R_buqiang+1
            s1 = m.ConstrainedSketch(name='R_strengtheningLayer'+str(R_buqiang), sheetSize=200.0)
            cline=s1.ConstructionLine((0,10),(0,-10))
            s1.assignCenterline(line=cline)

     
            #外轮廓在切线方向加上
            xvalue,yvalue=[],[]

            
            for i in range(len(fengtoui)-1):
                k1=(solvek1[i+1][1]-solvek1[i][1])/(solvek1[i+1][0]-solvek1[i][0])
                k2=-1/k1
                theta1=arctan(k2)#弧度值
                xvalue.append(fengtoui[i][0]+thick_list[i]*cos(theta1))
                yvalue.append(fengtoui[i][1]+thick_list[i]*sin(theta1))
            # xvalue.append(fengtoui[-1][0]+thick_list[-1]*cos(theta1))
            # yvalue.append(fengtoui[-1][0]+thick_list[-1]*sin(theta1))
            fengtouo=zip(xvalue,yvalue)
            
            for i in range(len(fengtoui)-1):
                s1.Line(point1=fengtoui[i],point2=fengtoui[i+1]) 

            for i in range(len(fengtouo)-1):
                s1.Line(point1=fengtouo[i],point2=fengtouo[i+1]) 
            
            #s1.Spline(points=fengtoui)    
            #s1.Spline(points=fengtouo)      
            
            #连线
            # line1=s1.Line(point1=(right_point[0],right_point[1]), point2=(fengtoui[-1][0],fengtoui[-1][1]))
            # line2=s1.Line(point1=(right_point[0],right_point[1]), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            line2=s1.Line(point1=(fengtoui[-1][0],fengtoui[-1][1]), point2=(fengtouo[-1][0],fengtouo[-1][1]))
            #line3=s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]), point2=(fengtouo[0][0],fengtouo[0][1]))
            line3_1 = s1.Line(point1=(fengtouo[0][0],fengtouo[0][1]),point2=(left_point[0],left_point[1]))
            line3_2 = s1.Line(point1=(fengtoui[0][0],fengtoui[0][1]),point2=(left_point[0],left_point[1]))

            
            R_lxcount=R_lxcount+1
            #更新参考轮廓
            if x2 >= L:
                cankaolunkuo1=kfengtoui[:indexnum_i]+fengtouo
                cankaolunkuo1.append(fengtoui[-1])
                cankaolunkuo=cankaolunkuo1+kfengtoui[indexnum_j:]
            else:
                cankaolunkuo=kfengtoui[:indexnum_i]+fengtouo+kfengtoui[indexnum_j:]
            
            p1 = m.Part(name='R_strengtheningLayer'+str(R_buqiang), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s1, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)

            #设置装配绑定参考面
           
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            luoxuanceng=a1.Instance(name='R_strengtheningLayer'+str(R_buqiang)+'-1', part=p1, dependent=ON)
            all_instances.append(luoxuanceng)
            all_name.append('R_strengtheningLayer'+str(R_buqiang)+'-1')

            #存储侧边参考面
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.pointOn
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            #设置材料赋予参考面
            f=p1.faces

            for i in range(len(fengtouo)-1):  
                side1Faces11 = f.findAt(((((fengtouo[i][0]+fengtouo[i+1][0])/2)*cos(0.5*rotangle*todegree),(fengtouo[i][1]+fengtouo[i+1][1])/2,-(fengtouo[i][0]+fengtouo[i+1][0])/2*sin(0.5*rotangle*todegree)), ))
                p1.Surface(side1Faces=side1Faces11, name='R_surfStrengtheningTop'+str(R_lxcount)+'-'+str(i+1))        
            

            
             
            p1 = m.parts['R_strengtheningLayer'+str(R_buqiang)]
            plant1=p1.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
            d1=p1.datums
            
            #封头画切分线
            t = p1.MakeSketchTransform(sketchPlane=f.findAt(coordinates=centercoordinate[0]
            ), sketchUpEdge=d1[plant1.id], sketchPlaneSide=SIDE1, origin=(0.0, 
                0.0, 0.0))
            s = m.ConstrainedSketch(name='__profile__', 
                sheetSize=1258.72, gridSpacing=31.46, transform=t)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=SUPERIMPOSE)
         
            p1.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
            
            pickface=f.findAt(coordinates=centercoordinate[0])
            
            
            for k in range(1,len(fengtoui)-2):
                s.Line(point1=fengtoui[k],point2=fengtouo[k])
            p1.PartitionFaceBySketch(sketchUpEdge=d1[plant1.id], faces=pickface,sketch=s)
            s.unsetPrimaryObject()
            del m.sketches['__profile__']            
            #封头处切分
            p1 = m.parts['R_strengtheningLayer'+str(R_buqiang)]
            c=p1.cells
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            for k in range(1,len(fengtoui)-2):
                c=p1.cells
                pickedEdges =(e1.findAt(coordinates=(0.5*(fengtoui[k][0]+fengtouo[k][0]), 0.5*(fengtoui[k][1]+fengtouo[k][1]), 0.0)),)   
                p1.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges)
            #设置集####且根据中心点坐标计算对应缠绕角度
            
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.getCentroid()
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0])
            face_coordinate.sort(key=(lambda x: x[0]))
            count=0
            list_angle=[]
            c=p1.cells
            if type_reinforce == "chanraobuqiang":
                for k in face_coordinate:
                    count=count+1
                    pickcells=c.findAt(((k[0], k[1], k[2]),))
                    #print r,r0,(k[0],k[1], k[2])
                    winding_angle=arcsin(r0/k[0])*180.0/pi
                    list_angle.append(winding_angle)
                    p1.Set(name='R_set-'+str(R_lxcount)+'-'+str(count),cells=pickcells)
                r_all_angle.append(list_angle)   
            else:
                for k in face_coordinate:
                    count=count+1
                    pickcells=c.findAt(((k[0], k[1], k[2]),))
                    winding_angle=buqiangangle
                    list_angle.append(winding_angle)
                    p1.Set(name='R_set-'+str(R_lxcount)+'-'+str(count),cells=pickcells)
                r_all_angle.append(list_angle)   
            #赋予角度
            r_all_count.append(count)


            #切分斜角前面多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[0], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[0][0], 
                fengtoui[0][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[0][0], fengtouo[0][1], 0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[0], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[0][0]+fengtouo[0][0]), 0.5*(fengtoui[0][1]+fengtouo[0][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], 
                -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges) 


    #        切分斜角后面多出来的一部分
            f1, e1, d2 = p1.faces, p1.edges, p1.datums
            pickedFaces = f1.findAt((face_coordinate[-1], ))
            v1, e1, d1 = p1.vertices, p1.edges, p1.datums
            p1.PartitionFaceByShortestPath(point1=v1.findAt(coordinates=(fengtoui[-2][0], 
                fengtoui[-2][1], 0.0)), point2=v1.findAt(coordinates=(fengtouo[-1][0], fengtouo[-1][1], 0.0)), faces=pickedFaces) 


            c = p1.cells
            pickedCells = c.findAt((face_coordinate[-1], ))
            e = p1.edges
            pickedEdges =(e.findAt(coordinates=(0.5*(fengtoui[-2][0]+fengtouo[-1][0]), 0.5*(fengtoui[-2][1]+fengtouo[-1][1]), 0.0)),) 
            p1.PartitionCellBySweepEdge(sweepPath=e.findAt(coordinates=(fengtouo[1][0]*cos(0.5*rotangle*todegree), fengtouo[1][1], 
                -1*fengtouo[1][0]*sin(0.5*rotangle*todegree))), cells=c, edges=pickedEdges) 
       

        

        else:           
            s2 = m.ConstrainedSketch(name='R_hoopLayer'+str(R_hxcount), sheetSize=200.0)
            cline=s2.ConstructionLine((0,10),(0,-10))
            #筒身
            tP = R_hoopThickness[R_hoopThicknessCount]
            tnewRi=R_thickness
            R_thickness+=tP
            tnewRo=tnewRi+tP
                      
            
            line1=s2.Line(point1=(tnewRi,y2), point2=(tnewRi,0))
            line2=s2.Line(point1=(tnewRo,y2), point2=(tnewRo,0))        
            #连线
            line1=s2.Line(point1=(cankaolunkuo[-1][0],cankaolunkuo[-1][1]), point2=(tnewRi,y2))
            line2=s2.Line(point1=(cankaolunkuo[-1][0],cankaolunkuo[-1][1]), point2=(tnewRo,y2))
            line3=s2.Line(point1=(tnewRi,0), point2=(tnewRo,0))
            R_hxcount=R_hxcount+1
            p1 = m.Part(name='R_hoopLayer'+str(R_hxcount), dimensionality=THREE_D, type=DEFORMABLE_BODY)
            p1.BaseSolidRevolve(sketch=s2, angle=rotangle, flipRevolveDirection=ON)
            #p1.setValues(geometryRefinement=EXTRA_FINE)
            part_name.append(p1)
            #设置材料赋予参考面
            f=p1.faces
            side1Faces1 = f.findAt(((tnewRi*cos(0.5*rotangle*todegree),y2/2.0,-tnewRi*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces1, name='R_surfHoopBottom')       
            side1Faces2 = f.findAt(((tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree)), ))
            p1.Surface(side1Faces=side1Faces2, name='R_surfHoopTop'+str(R_hxcount))       
            #设置装配绑定参考面 
            a1 = m.rootAssembly
            a1.DatumCsysByDefault(CARTESIAN)
            huanxiangceng=a1.Instance(name='R_hoopLayer'+str(R_hxcount)+'-1', part=p1, dependent=ON)  
            all_instances.append(huanxiangceng)
            all_name.append('R_hoopLayer'+str(R_hxcount)+'-1')
            c=p1.cells
            f_all=p1.faces
            face_coordinate=[]
            for face in f_all:
                centercoordinate=face.pointOn
                if centercoordinate[0][2] ==0.0:
                    face_coordinate.append(centercoordinate[0]) 
            for k in face_coordinate:
                pickcells=c.findAt(((k[0], k[1], k[2]),))
                p1.Set(name='R_set-'+str(R_hxcount),cells=pickcells)
             

            SDcoordinates=(tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree))
            # c = p1.cells
            # f = p1.faces
            # p1.assignStackDirection(referenceRegion=f.findAt(coordinates=(tnewRo*cos(0.5*rotangle*todegree),y2/2.0,-tnewRo*sin(0.5*rotangle*todegree))), cells=c)






    a1 = m.rootAssembly
    if len(all_instances)> 1:
        a1.InstanceFromBooleanMerge(name='Part-F', instances=all_instances, keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)
        p2 = m.parts['Part-F']
    else:
        msg = "set two layers at least" 
        raise AbaqusException, msg


###########删除组合前的部件
    for  i in range(L_lxcount-L_buqiang):
        del m.parts['L_helixLayer'+str(i+1)]
        del m.sketches['L_helixLayer'+str(i)]
    for i in range(L_hxcount): 
        del m.parts['L_hoopLayer'+str(i+1)]   
        del m.sketches['L_hoopLayer'+str(i)] 
    for i in range(L_buqiang): 
        del m.parts['L_strengtheningLayer'+str(i+1)]
        del m.sketches['L_strengtheningLayer'+str(i+1)]
        
    for  i in range(R_lxcount-R_buqiang):
        del m.parts['R_helixLayer'+str(i+1)]
        del m.sketches['R_helixLayer'+str(i)]
    for i in range(R_hxcount): 
        del m.parts['R_hoopLayer'+str(i+1)]
        del m.sketches['R_hoopLayer'+str(i)] 
    for i in range(R_buqiang): 
        del m.parts['R_strengtheningLayer'+str(i+1)]
        del m.sketches['R_strengtheningLayer'+str(i+1)]


    #del a1.features['Part-1-1']
    


    ##########重新赋予角度##########
    
    L_lxcount=0
    L_hxcount=0
    L_buqiang=0
    for windtype in L_tables:
        if windtype == 'helix':
            L_lxcount=L_lxcount+1
            scount=all_count[L_lxcount-1]

            list_angle=all_angle[L_lxcount-1]
            
            for k in range(1,scount):
                layupOrientation = None
                region1 = p2.sets['L_set-'+str(L_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['L_surfHeadTop'+str(L_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='L_helix-'+str(L_lxcount)+'-head-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 
            

            
            
            
            #筒身段赋予角度
            layupOrientation = None
            region1 = p2.sets['L_set-'+str(L_lxcount)+'-'+str(scount)]
            normalAxisRegion = p2.surfaces['L_surfCylinderTop'+str(L_lxcount)]
            compositeLayup = p2.CompositeLayup(
                name='L_helix-'+str(L_lxcount)+'-cylinder', description='', elementType=SOLID, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)    

        elif windtype == 'helixReaming':  
            L_lxcount=L_lxcount+1
            scount=all_count[L_lxcount-1]
            list_angle=all_angle[L_lxcount-1]    
            for k in range(1,scount):
                layupOrientation = None
                region1 = p2.sets['L_set-'+str(L_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['L_surfHeadTop'+str(L_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='L_helix-'+str(L_lxcount)+'-head-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 

            #筒身段赋予角度
            layupOrientation = None
            region1 = p2.sets['L_set-'+str(L_lxcount)+'-'+str(scount)]
            normalAxisRegion = p2.surfaces['L_surfCylinderTop'+str(L_lxcount)]
            compositeLayup = p2.CompositeLayup(
                name='L_helix-'+str(L_lxcount)+'-cylinder', description='', elementType=CONTINUUM_SHELL, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)  

            
        elif windtype == "strengthening":  
            L_buqiang=L_buqiang+1
            L_lxcount=L_lxcount+1
            scount=all_count[L_lxcount-1]
            list_angle=all_angle[L_lxcount-1]  
           
            for k in range(1,scount+1):
                layupOrientation = None
                region1 = p2.sets['L_set-'+str(L_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['L_surfStrengtheningTop'+str(L_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='L_strengthening-'+str(L_buqiang)+'-strengthening-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 
        else:    
            #环向层赋予角度
            L_hxcount=L_hxcount+1
            layupOrientation = None
            region1 = p2.sets['L_set-'+str(L_hxcount)]
            normalAxisRegion = p2.surfaces['L_surfHoopTop'+str(L_hxcount)]
            compositeLayup = p2.CompositeLayup(
                name='L_hoop-'+str(L_hxcount), description='', elementType=SOLID, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=huanxiang_alpha, 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            # compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                # material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                # orientationType=SPECIFY_ORIENT, orientationValue=-1*huanxiang_alpha, 
                # additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                # axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)

    ##########################################################################################################
    ##########################################################################################################

    R_lxcount=0
    R_hxcount=0
    R_buqiang=0
    for windtype in R_tables:
        if windtype == 'helix':
            R_lxcount=R_lxcount+1
            scount=r_all_count[R_lxcount-1]

            list_angle=r_all_angle[R_lxcount-1]
            for k in range(1,scount):
                layupOrientation = None
                region1 = p2.sets['R_set-'+str(R_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['R_surfHeadTop'+str(R_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='R_helix-'+str(R_lxcount)+'-head-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 
            

            
            
            
            #筒身段赋予角度
            layupOrientation = None
            region1 = p2.sets['R_set-'+str(R_lxcount)+'-'+str(scount)]
            normalAxisRegion = p2.surfaces['R_surfCylinderTop'+str(R_lxcount)]
            compositeLayup = p2.CompositeLayup(
                name='R_helix-'+str(R_lxcount)+'-cylinder', description='', elementType=SOLID, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)    

        elif windtype == 'helixReaming':  
            R_lxcount=R_lxcount+1
            scount=r_all_count[R_lxcount-1]
            list_angle=r_all_angle[R_lxcount-1]    
            for k in range(1,scount):
                layupOrientation = None
                region1 = p2.sets['R_set-'+str(R_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['R_surfHeadTop'+str(R_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='R_helix-'+str(R_lxcount)+'-head-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 

            #筒身段赋予角度
            layupOrientation = None
            region1 = p2.sets['R_set-'+str(R_lxcount)+'-'+str(scount)]
            normalAxisRegion = p2.surfaces['R_surfCylinderTop'+str(R_lxcount)]
            compositeLayup = p2.CompositeLayup(
                name='R_helix-'+str(R_lxcount)+'-cylinder', description='', elementType=CONTINUUM_SHELL, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[-1], 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)  

            
        elif windtype == "strengthening":  
            R_buqiang=R_buqiang+1
            R_lxcount=R_lxcount+1
            scount=r_all_count[R_lxcount-1]
            list_angle=r_all_angle[R_lxcount-1]   
            for k in range(1,scount+1):
                layupOrientation = None
                region1 = p2.sets['R_set-'+str(R_lxcount)+'-'+str(k)]
                normalAxisRegion = p2.surfaces['R_surfStrengtheningTop'+str(R_lxcount)+'-'+str(k)]
                compositeLayup = p2.CompositeLayup(
                    name='R_strengthening-'+str(R_buqiang)+'-strengthening-'+str(k), description='', elementType=SOLID, 
                    symmetric=False)
                compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                    useDensity=OFF)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                    material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                    additionalRotationType=ROTATION_NONE, angle=0.0, 
                    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                    normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                    normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                    primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                    primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 
        else:    
            #环向层赋予角度
            R_hxcount=R_hxcount+1
            layupOrientation = None
            region1 = p2.sets['R_set-'+str(R_hxcount)]
            normalAxisRegion = p2.surfaces['R_surfHoopTop'+str(R_hxcount)]
            compositeLayup = p2.CompositeLayup(
                name='R_hoop-'+str(R_hxcount), description='', elementType=SOLID, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
                material=compositeName, thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                orientationType=SPECIFY_ORIENT, orientationValue=huanxiang_alpha, 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)
            # compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region1, 
                # material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                # orientationType=SPECIFY_ORIENT, orientationValue=-1*huanxiang_alpha, 
                # additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                # axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(0.0, 1.0, 0.0), 
                primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)





    #控制网格
    ii=p2.vertices
    e1=p2.edges
    meshEdgeCount = 10
    for edge in e1:
        zu=edge.getVertices()
        zu1=ii[zu[0]].pointOn
        zu2=ii[zu[1]].pointOn
        p1x,p1y,p1z,p2x,p2y,p2z=zu1[0][0],zu1[0][1],zu1[0][2],zu2[0][0],zu2[0][1],zu2[0][2]
        distance=((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2)**0.2
        if zu1[0][1]==zu2[0][1] and zu1[0][2]*zu2[0][2]==0 and zu1[0][2]+zu2[0][2]!=0:
            meshEdgeCount-=1
            sd=edge.pointOn
            pickedges=e1.findAt(sd)
            p2.seedEdgeByNumber(edges=pickedges, number=piece_number, constraint=FINER)
            if meshEdgeCount==0:
                break
     
    p2.seedPart(size=approximate_size, deviationFactor=0.1, minSizeFactor=0.1)



    c = p2.cells
    f = p2.faces
    p2.assignStackDirection(referenceRegion=f.findAt(coordinates=SDcoordinates), cells=c)

    upFlag=False
    downFlag=False
    
    ds=p2.datums
    for x in L_style:
        if x[5] == "up":
            upFlag=True
        if x[5] == "down":
            downFlag=True
    
    if upFlag:
        PL1=p2.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=leftAcc)
        c = p2.cells
        p2.PartitionCellByDatumPlane(cells=c,datumPlane=ds[PL1.id])
    if downFlag:
        PL2=p2.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=rightAcc)
        c = p2.cells
        p2.PartitionCellByDatumPlane(cells=c,datumPlane=ds[PL2.id])





    p2.generateMesh()

    a1.ReferencePoint(point=(0.0, 0.0, 0.0))
    a1.ReferencePoint(point=(0.0, 200, 0.0))


    m.StaticStep(name='Step-1', previous='Initial')
    m.steps['Step-1'].setValues(nlgeom=ON)     ###几何非线性
    m.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U'))


    field_count=1



    #########对称边界条件


            
            


    """

    sub=[1]
    if sub[0]==0:
        E1=114000
        E2=8610
        E3=8610
        V12=0.3
        V13=0.3
        V23=0.45
        G12=4160
        G13=4160
        G23=3000

        mat=[E1,E2,E3,V12,V13,V23,G12,G13,G23] #初始材料参数
        mt=[1*E1,    0.2*E2,  1*E3,    1*V12,    1*V13,    1*V23,    0.2*G12,  1*G13,    0.2*G23  ] #基体拉伸
        mc=[1*E1,    0.4*E2,  1*E3,    1*V12,    1*V13,    1*V23,    0.4*G12,  1*G13,    0.4*G23  ] #基体压缩
        ft=[0.07*E1, 0.07*E2, 0.07*E3, 0.07*V12, 0.07*V13, 0.07*V23, 0.07*G12, 0.07*G13, 0.07*G23 ] #纤维拉伸
        fc=[0.14*E1, 0.14*E2, 0.14*E3, 0.14*V12, 0.14*V13, 0.14*V23, 0.14*G12, 0.14*G13, 0.14*G23 ] #纤维压缩
        gs=[1*E1,    1*E2,    1*E3,    0.001*V12,    1*V13,    1*V23,    0.001*G12,    1*G13,    1*G23    ] #剪切
        ff=[1*E1,    1*E2,    0.001*E3,    1*V12,    0.001*V13,    0.001*V23,    1*G12,    0.001*G13,    0.001*G23    ] #分层

        dt=[ft,fc,mt,mc,gs,ff]




        m=[0,0,0,0,0,0]
        constants=[]
        #No damage
        d=[E1,E2,E3,V12,V13,V23,G12,G13,G23, 0, 0, 0, 0, 0, 0]
        #print str(tuple(d))
        a0=(tuple(d))
        constants.append(a0)



        #one damge
        i=0
        while(i<6):
            d=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for k in range(9):
                d[k]=round(dt[i][k],4)
                d[9+i]=1
            #print str(tuple(d))
            a1=(tuple(d))
            constants.append(a1)
            i+=1

        #two damge
        i=0
        while(i<6):    
            j=i+1
            while(j<6):   
                d=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                for k in range(9):
                    d[k]=round(min(dt[i][k],dt[j][k]),4)
                    d[9+i]=1
                    d[9+j]=1
                #print str(tuple(d))
                a2=(tuple(d))
                constants.append(a2)        
                j+=1
            i+=1
        #three damage
        i=0
        while(i<6):    
            j=i+1
            while(j<6):
                l=j+1
                while(l<6):
                    d=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                    for k in range(9):
                        d[k]=round(min(dt[i][k],dt[j][k],dt[l][k]),4)
                        d[9+i]=1
                        d[9+j]=1
                        d[9+l]=1
                    #print str(tuple(d))
                    a3=(tuple(d))
                    constants.append(a3)
                    l+=1
                j+=1
            i+=1
        #for damage
        i=0
        while(i<6):    
            j=i+1
            while(j<6):
                l=j+1
                while(l<6):
                    m=l+1
                    while(m<6):
                        d=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                        for k in range(9):
                            d[k]=round(min(dt[i][k],dt[j][k],dt[l][k],dt[l][k]),4)
                            d[9+i]=1
                            d[9+j]=1
                            d[9+l]=1
                            d[9+m]=1
                        #print str(tuple(d))
                        a4=(tuple(d))
                        constants.append(a4)
                        m+=1
                    l+=1
                j+=1
            i+=1    
        #five damage
        i=0
        while(i<6):    
            j=i+1
            while(j<6):
                l=j+1
                while(l<6):
                    m=l+1
                    while(m<6):
                        p=m+1
                        while(p<6):
                            d=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                            for k in range(9):
                                d[k]=round(min(dt[i][k],dt[j][k],dt[l][k],dt[l][k],dt[p][k]),4)
                                d[9+i]=1
                                d[9+j]=1
                                d[9+l]=1
                                d[9+m]=1
                                d[9+p]=1
                            #print str(tuple(d)) 
                            a5=(tuple(d)) 
                            constants.append(a5)
                            p+=1
                        m+=1
                    l+=1
                j+=1
            i+=1
        #six damage
        d=[0.07*E1, 0.07*E2, 0.07*E3, 0.07*V12, 0.07*V13, 0.07*V23, 0.07*G12, 0.07*G13, 0.07*G23 , 1, 1, 1, 1, 1, 1]
        #print str(tuple(d)) 
        a6=(tuple(d))
        constants.append(a6)




        mat_comp.Depvar(n=6)
        mat_comp.UserDefinedField()
        mat_comp.Elastic(dependencies=6, table=(constants), type=ENGINEERING_CONSTANTS)
    """
    
  
    
    
    
    
    
    
    
    
    
    
    
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 


# 下面是压力容器的
######################################################################################################################

# 自定义数据处理相关函数
#声明自定义函数
def ChaiFen(data):
    a = data
    Line = zeros((1, 4))
    Aa = zeros((len(a), 4))
    for i in range(len(a)):  # 循环赋值
        Line0 = a[i].split()
        Aa[i][0] = float(Line0[0])
        Aa[i][1] = float(Line0[1])/2#将文件中的直径值转换为半径值
        Aa[i][2] = float(Line0[2])
        Aa[i][3] = float(Line0[3])
    Aa = Aa[Aa[:, 0].argsort()]  # 按照第1列对行排序
       # Rmax = max(Aa[0:, 1])  # 获取最大值
    # tmin = min(Aa[0:, 3])  # 获取最小值
    # c = np.where(Aa == max(Aa[0:, 1]))  # 获取直径最大值所在行号和列号，即是筒身段的直径值
    # return [Aa,c,Rmax,tmin]
    return [Aa]

def DataSolve(data):
    chaifen0=ChaiFen(data)#拆分出各个值
    Aa=chaifen0[0]

#预定义
    # 轮廓绘制用
    X =[]
    R =[]
    alpha =[]
    h = []
    # 角度和厚度用
    X11 = []
    R11 = []
    alpha11 = []
    h11 = []

#循环赋值
    for i1 in range(len(Aa)):
        # 画轮廓
        X.append(Aa[i1][0]);
        R.append(Aa[i1][1]);#半径值，文件里面是直径值，已在拆分函数里面转换为半径
        # alpha.append(Aa[i1][2]);
        # h.append(Aa[i1][3]);

        # 赋角度和厚度
        # if Aa[i1][2]>0:#只提取一个循环的数据
        #     X11.append(Aa[i1][0]);
        #     R11.append(Aa[i1][1]);  # 半径值，文件里面是直径值，已在拆分函数里面转换为半径
        #     alpha11.append(Aa[i1][2]);
        #     h11.append(Aa[i1][3]);

        # # 输出全部厚度
        X11.append(Aa[i1][0]);
        R11.append(Aa[i1][1]);  # 半径值，文件里面是直径值，已在拆分函数里面转换为半径
        alpha11.append(Aa[i1][2]);
        h11.append(Aa[i1][3]);

# 获取重复值行号
    hanghao= []
    for ii1 in range(len(Aa) - 1):
        if Aa[ii1 + 1][0] == Aa[ii1][0] and abs(Aa[ii1 + 1][2]) == abs(Aa[ii1][2]):  # 位同角同(取绝对值)
            hanghao.append(ii1)
        else:
            pass
        if Aa[ii1 + 1][0] == Aa[ii1][0] and abs(Aa[ii1 + 1][2]) != abs(Aa[ii1][2]):  # 位同角不同(取绝对值)
            hanghao.append(ii1)
        else:
            pass

    hanghao1 = []
    for ii1 in range(len(X11) - 1):
        if X11[ii1 + 1] == X11[ii1] and abs(alpha11[ii1 + 1]) == abs(alpha11[ii1]):  # 位同角同(取绝对值)
            hanghao1.append(ii1)
        else:
            pass
        if X11[ii1 + 1] == X11[ii1] and abs(alpha11[ii1 + 1]) != abs(alpha11[ii1]):  # 位同角不同(取绝对值)
            hanghao1.append(ii1)
        else:
            pass

# 去除重复坐标值
    i=0
    for num in hanghao:
        num=num-i
        i=i+1
        del (X[num])
        del (R[num])
    i = 0
    for num in hanghao1:
        num = num - i
        i = i + 1
        del (X11[num])
        del (R11[num])
        del (alpha11[num])
        del (h11[num])


    return [X,R,X11,R11,alpha11,h11]###X11,alpha11原始LAM文件位置角度

 #BSpline3次
# degree. 3次样条
    # 第一步：设置曲线次数为p = 3，读取控制顶点P0, P1, ...Pn，根据控制顶点得到节点向量U = {u0, u1....um}，注意m = n + p + 1.
    # 为了保证拟合的曲线经过第一个、最后一个控制点，节点向量首末重复度要设置为为p + 1，即U = {0, 0, 0, 0, u(4)...u(n), 1, 1, 1, 1}
    # 第二步：按照本文描述的DeBoor算法完成递推函数
    # 第三步：在0 ≤ u ≤ 1取值，即可得到u对应的P(u)位置
# 规范化初始弧长
def guifanhua(CPointX, CPointY):
    knotType = 1  # knot的计算方式 0：平均取值  1：根据长度取值
    k = 3
    n = len(CPointX)
    knot = list(nmp.arange(n + k + 1))  # 创建首末控制点重度，便于曲线过端点控制点
    for i in range(0, 4):  # 由B样条性质构建左右数据控制点的k+1重节点，k为基函数阶数
        knot[i] = 0.0
        knot[n + i] = 1.0
    if knotType:
        L = list(nmp.arange(n - 1))
        S = 0
        for i in range(n - 1):  # 计算矢量长
            L[i] = nmp.sqrt(pow(CPointX[i + 1] - CPointX[i], 2) + pow(CPointY[i + 1] - CPointY[i], 2))
            S = S + L[i]
        tmp = L[0]
        for i in range(4, n):  # 按矢量比例均分，曲线定义域内
            tmp = tmp + L[i - 3]
            knot[i] = tmp / S
    else:
        tmp = 1 / (float(n) - 3)  # 因为三次样条是4个节点一个曲线段，现在在将曲线定义域用n-k来均分定义区间
        tmpS = tmp
        for i in range(4, n):
            knot[i] = tmpS
            tmpS = tmpS + tmp
    return knot
# 基函数值递归计算
def coxDeBoor(u, knots, i, k):
    # Test for end conditions
    if (k == 0):
        if (knots[i] <= u and u < knots[i + 1]):
            return 1
        return 0
    Den1 = knots[i + k] - knots[i]
    Den2 = knots[i + k + 1] - knots[i + 1]
    Eq1 = 0;
    Eq2 = 0;
    if Den1 > 0:
        Eq1 = ((u - knots[i]) / Den1) * coxDeBoor(u, knots, i, (k - 1))
    if Den2 > 0:
        Eq2 = ((knots[i + k + 1] - u) / Den2) * coxDeBoor(u, knots, (i + 1), (k - 1))

    return Eq1 + Eq2
# 生成Bspline点坐标，并修正R值
def B3Spline(CPointX,CPointY,knot,R0_jikong2,Bsjindu):
    BsplineX = []
    BsplineY = []
    for ip in range(len(CPointX)-3):#创建B样条
        u =linspace(knot[ip+3], knot[ip + 4],Bsjindu)#节点区间值
        for j in arange(len(u)):#创建均B样条
            BsplineX.append(CPointX[ip]*coxDeBoor(u[j], knot, ip, 3)+CPointX[ip+1]*coxDeBoor(u[j], knot, ip+1, 3)\
                +CPointX[ip+2]*coxDeBoor(u[j], knot, ip+2, 3)+CPointX[ip+3]*coxDeBoor(u[j], knot, ip+3, 3))

            BsplineY.append(CPointY[ip] * coxDeBoor(u[j], knot, ip, 3) + CPointY[ip + 1] * coxDeBoor(u[j], knot, ip + 1, 3) \
                 + CPointY[ip + 2] * coxDeBoor(u[j], knot, ip + 2, 3) + CPointY[ip + 3] * coxDeBoor(u[j], knot, ip + 3, 3))
    # 2020.8.30改
    # 会导致出现点过原点的错误
    if BsplineX[-1]==0:
        del BsplineX[-1]
        del BsplineY[-1]

    x1 = BsplineX[-1] + (R0_jikong2 - BsplineY[-2]) * ((BsplineX[-1] - BsplineX[-2]) / (BsplineY[-1] - BsplineY[-2]))
    BsplineX.append(x1)
    BsplineY.append(R0_jikong2)
    # plt.plot(BsplineX, BsplineY, '-r')
    # plt.plot([BsplineX[0], BsplineX[-1]], [0, 0], 'g')
    # plt.plot([BsplineX[0], BsplineX[-1]], [R0_jikong, R0_jikong], 'b')
    # plt.plot(CPointX, CPointY, '-ob', label="Waypoints")
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()
    return [BsplineX,BsplineY]

#通过厚度求取外轮廓，即是通过两个数据点的向量值夹角和厚度值反求对应的外轮廓点
def wailunkuo1(X,R,H):
    X2_zuo_wai = []
    R2_zuo_wai = []

# 第一个数据点位的厚度反算轮廓
    theta = math.atan(-(X[1] - X[0]) / (R[1] - R[0]))
    Delta_X = (2*H[0]) * cos(theta)
    Delta_R = (2*H[0]) * sin(theta)
    X2_zuo_wai.append(X[0] - Delta_X)
    R2_zuo_wai.append(R[0] - Delta_R)

# 除第一个数据点外的数据点位厚度反算轮廓
    for i111 in range(len(X) - 1):
        i111 = i111 + 1
        theta=math.atan(-(R[i111]-R[i111-1])/(X[i111]-X[i111-1]))
        Delta_X=(2*H[i111])*sin(theta)
        Delta_R=(2*H[i111])*cos(theta)
        if R[i111]>R[i111-1]:
            X2_zuo_wai.append(X[i111]+Delta_X)
            R2_zuo_wai.append(R[i111]+Delta_R)
        if R[i111]<R[i111-1]:
            X2_zuo_wai.append(X[i111]+Delta_X)
            R2_zuo_wai.append(R[i111]+Delta_R)
        if R[i111]==R[i111-1]:
            X2_zuo_wai.append(X[i111])
            R2_zuo_wai.append(R[i111]+2*H[i111])

    # plt.plot(X2_zuo_wai,R2_zuo_wai,'-ob', label="Waypoints")
    # plt.show()
    return [X2_zuo_wai,R2_zuo_wai]

#通过厚度求取外轮廓，即是通过两个数据点的向量值夹角和厚度值反求对应的外轮廓点
#左封头
def wailunkuozuo(X,R,H,tsX):
    X2_zuo_wai = list()
    R2_zuo_wai = list()
    # if X[0]!=max(X):
    #     X=X[::-1]#因为以赤道圆作为初始位置，需翻转左封头的数据
    #     R=R[::-1]
    #     H=H[::-1]
    X2_zuo_wai.append(tsX[0])
    R2_zuo_wai.append(max(R)+2*H[0])
    for i111 in range(len(X) - 1):
        i111 = i111 + 1
        theta=math.atan(-(X[i111]-X[i111-1])/(R[i111]-R[i111-1]))
        Delta_X=(2*H[i111])*cos(theta)
        Delta_R=(2*H[i111])*sin(theta)
        X2_zuo_wai.append(X[i111]-Delta_X)
        R2_zuo_wai.append(R[i111]-Delta_R)
    # plt.plot(X2_zuo_wai,R2_zuo_wai)
    return [X2_zuo_wai,R2_zuo_wai]

#右封头
def wailunkuoyou(X,R,H,tsX):
    X2_you_wai = list()
    R2_you_wai = list()

    X2_you_wai.append(tsX[1])
    R2_you_wai.append(max(R)+ 2*H[0])
    for i222 in range(len(X) - 1):
        i222 = i222 + 1
        theta=math.atan(-(X[i222]-X[i222-1])/(R[i222]-R[i222-1]))
        Delta_X=(2*H[i222])*cos(theta)
        Delta_R=(2*H[i222])*sin(theta)
        X2_you_wai.append(X[i222]+Delta_X)
        R2_you_wai.append(R[i222]+Delta_R)
    # plt.plot(X2_you_wai,R2_you_wai)
    return [X2_you_wai,R2_you_wai]

# 原始数据和B样条之间的修正替代
def XiuZheng(x,X,R):#R坐标修正
    # 修正半径值,因为X/Y坐标是一一对应的，所以可以找到离x坐标最近的那个R作为我的现在的修正的R,为后面厚度反推外轮廓做准备
    Rxiuzh = []
    Xxiuzh = []
    Hxiuzh = []
    for i in range(len(x)):
        jdz = [abs(X[j] - x[i]) for j in range(len(X))]  # 绝对值
        c1 = jdz.index(min(jdz))  # 获取各个x对应为修正位置
        Xxiuzh.append(X[c1])
        Rxiuzh.append(R[c1])
    return [Xxiuzh,Rxiuzh]

def houdupingwen(X,R,pingwenxishu,tsX):#前后前跳动值平稳处理，用来为轮廓生成做准备
    # 两个向量夹角判定
    j=0
    bb=[]
    weizhi=[]
    for i in range(len(X) - 2):
        axaingliang=[X[i+1-j]-X[i-j],R[i+1-j]-R[i-j]]#a向量
        bxaingliang=[[X[i+2-j]-X[i+1-j]],[R[i+2-j]-R[i+1-j]]]#b向量
        amochang=((X[i+1-j]-X[i-j])**2+(R[i+1-j]-R[i-j])**2)**0.5#模长
        bmochang=((X[i+2-j]-X[i+1-j])**2+(R[i+2-j]-R[i+1-j])**2)**0.5#模长
        alp=(math.acos((np.dot(axaingliang,bxaingliang)) / (amochang* bmochang)))*180/math.pi
        bb.append(alp)
        if abs(alp)>pingwenxishu*180/math.pi and X[i-j]<tsX[0]:#去除急转的角度相邻位置值
            weizhi.append(i+2)#将波动位置存储起来后面统一删除
        if abs(alp)>pingwenxishu*180/math.pi and X[i-j]>tsX[1]:#去除急转的角度相邻位置值
            weizhi.append(i+2)#将波动位置存储起来后面统一删除
    for ii in range(len(weizhi)):
        del (X[weizhi[ii]-j])
        del (R[weizhi[ii]-j])
        j+=1
    # for ij in arange(0,10,1):#平稳迭代次数
    #     alp = []
    #     for i in range(len(X)-1):
    #         alp.append(np.arctan((R[i+1]-R[i])/(X[i+1]-X[i])))
    #     weizhi=[i for i in range(len(alp)-1) if abs(alp[i+1])-abs(alp[i])>pingwenxishu]#去除急转的角度相邻位置值
    #     for i in range(len(weizhi)):
    #         del (X[weizhi[i]+1-i])
    #         del (R[weizhi[i]+1-i])
    return [X,R]

def Qujiange(old_listX,old_listR,shuliang):
    new_listX=[]
    new_listR=[]
    new_listX.append(old_listX[0])
    new_listR.append(old_listR[0])
    # 获取平均距离
    distance=[((old_listX[i]-old_listX[i+1])**2+(old_listR[i]-old_listR[i+1])**2)**0.5 for i in range(len(old_listX)-1)]
    pingjudistance=sum(distance)/shuliang
    i = 0
    for j in range(1,len(old_listX)):
        if ((old_listX[j]-old_listX[i])**2+(old_listR[j]-old_listR[i])**2)**0.5>pingjudistance:
            new_listX.append(old_listX[j-1])
            new_listR.append(old_listR[j-1])
            i=j
    # new_list=old_list[0:len(old_list):inter]
    if new_listX[-1]!=old_listX[-1] or new_listR[-1]!=old_listR[-1]:
        new_listX.append(old_listX[-1])
        new_listR.append(old_listR[-1])
    else:
        pass
    return [new_listX,new_listR]

# 2020.10.30
# 通过内部层
# 2020.8.23修正增加，这里的x会出现相等的情况必须再次去重处理
def quchong(X,R,H,alpha1):
    weizhi=[i for i in range(len(X)) if X[0:i+1].count(X[i])>1]# 找到去除X重复项的位置
    for i in range(len(weizhi)):
        del (X[weizhi[i]-i])
        del (R[weizhi[i]-i])
        del (H[weizhi[i]-i])
        del (alpha1[weizhi[i]-i])
    return [X,R,H,alpha1]

# 厚度光滑性处理
# def smooth(H):
#     # 对厚度进行平稳处理
#     h12=[]
#     h13=[]
#     # 左
#     for i in range(len(H)):
#         if X12[i]<0:
#             h12.append(sum(H[i:i+5])/5)
#     # 右
#     X121=list(reversed(X12))
#     h11=list(reversed(H))
#     for i in range(len(H)):
#         if X121[i] > tschangdu:
#             h13.append(sum(h11[i:i+5]) / 5)
#     # 重新拼装厚度列表
#     # 左
#     for i in range(len(h12)):
#         H[i]=h12[i]
#     # 右
#     for i in range(len(h13)):
#         H[-(i+1)]=h13[i]
#     return H



######################################################################################################################
 
 ######################################################################################################################
# def createModel(path,layUp,E1,E2,E3,v12,v13,v23,G12,G13,G23,pressureValue,hoopSingleThick,
        # L_polorHole,R_polarHole,totalLength,L_headLength,R_headLength,cylinderLength,rotAngle):
def createModel(path,layUp,E1,E2,E3,v12,v13,v23,G12,G13,G23,pressureValue,hoopSingleThick,
        cylinderLength,rotAngle):
    print(layUp)
    print(E1,E2,E3,v12,v13,v23,G12,G13,G23,pressureValue,hoopSingleThick,
        cylinderLength,rotAngle)
    
    # 导入相关的cadwind数据
    # 读取数据文件
    Lam1 = open(path, 'r')
    data1 = Lam1.readlines()
    Lam1.close()
    del (data1[0:2])  # 删除前面的两行非数值项,Cadwind默认生成的非数据行


    ######################################################################################################################


    puceng=[]

    for x in layUp:
        puceng.append(x[0])

        


    #自定义铺层顺序
    #puceng=['lx','hx','lx','hx','lx','hx']
    engineeringdata=(E1,E2,E3,v12,v13,v23,G12,G13,G23) 
    pressure=pressureValue#设计爆压
    hxdt=hoopSingleThick#环向单层厚度
    # lxdt=helixSingleThick#螺旋单层厚度
    # luoxuanjiao=helixAngle
    # R0_jikong1=L_polorHole#定义左侧极孔半径，便于求交点
    # R0_jikong2=R_polarHole#定义右侧极孔半径，便于求交点
    # changdu=totalLength#总长度
    # fengtouchazuo=L_headLength#左封头长度
    # fengtouchayou=R_headLength#右封头长度
    tschangdu=cylinderLength#筒身长度
    zhuanjiao=rotAngle#三维实体模型的转角


    #默认参数不需要插件输入
    #Bsjindu=BsplinePrecision#B样条单个区间拟合精度   
    # todegree=pi/180.0
    # jiangeshu=internalNum#绘制三维模型间隔取点数量
    # pingwenxishu=smoothCoefficient#外轮廓相邻点平稳系数，越小越平稳，太小会出问题，需适当调节

    Bsjindu=200#B样条单个区间拟合精度
    zhuanjiao=30.0#三维实体模型的转角
    todegree=pi/180.0
    pingwenxishu=0.2#外轮廓相邻点平稳系数，越小越平稳，太小会出问题，需适当调节


    ######################################################################################################################


    # 初始化数据存储参数
    # 端点值
    # 左
    XddianZ=[]
    RddianZ=[]
    # 右
    XddianY=[]
    RddianY=[]

    # 轮廓
    # 内表面
    neisX=[]
    neisR=[]
    # 外表面
    waisX=[]
    waisR=[]



       
     # 得到传入处理后的初步数据
    #内表面
    # 得到最内层的轮廓和螺旋缠绕的角度数据
    # 得到各个数据值
    data=DataSolve(data1)
    #绘制轮廓用
    X1=data[0]
    R1=data[1]
    # 角度赋予用
    X12=data[2]
    R12=data[3]
    alpha1=data[4]
    h1=data[5]


######################################################################################################################


    # 2020.11.1
    # 极孔线性外插
    # 左
    x1=X12[0]
    x2=X12[1]

    y1=R12[0]
    y2=R12[1]

    Y=np.array([[y1],[y2]])
    mat2=np.array([[x1,1],[x2,1]])
    mat2 = np.mat(mat2)
    mat3 = mat2.I#求逆
    xishu=mat3*Y
    zuojikquzhengX=-ceil(abs(X12[0]))
    zuojikquzhengR=round(xishu[0,0]*zuojikquzhengX+xishu[1,0])

    R0_jikong1=zuojikquzhengR
    fengtouchazuo=abs(zuojikquzhengX)#左封头长度

    # 右
    x1=X12[-1]
    x2=X12[-2]

    y1=R12[-1]
    y2=R12[-2]

    Y=np.array([[y1],[y2]])
    mat2=np.array([[x1,1],[x2,1]])
    mat2 = np.mat(mat2)
    mat3 = mat2.I
    xishu=mat3*Y
    youjikquzhengX=ceil(abs(X12[-1]))
    youjikquzhengR=round(xishu[0,0]*youjikquzhengX+xishu[1,0])

    R0_jikong2=youjikquzhengR
    fengtouchayou=youjikquzhengX-tschangdu#右封头长度

    # 2020.11.1
    # 找出筒身长度
    # Xts=[]
    # for i in range(len(X1)-2):
    #     if R1[i]==R1[i+1]==R1[i+2] or R1[i]==R1[i-1]==R1[i-2]:
    #         Xts.append(X1[i])
    # tschangdu=Xts[-1]-Xts[0]


    # 根据两边极孔条件删除不在范围内的点
    weizhi00=[i for i in range(len(R12)) if R12[i]<R0_jikong1 and X12[i]<0]
    j=0
    for i in range(len(weizhi00)):
        del (X12[weizhi00[i-j]])
        del (R12[weizhi00[i-j]])
        del (alpha1[weizhi00[i-j]])
        del (h1[weizhi00[i-j]])
        j+=1
        # weizhi00=[weizhi00[i]-1 for i in range(len(weizhi00))]
    weizhi00=[i for i in range(len(R12)) if R12[i]<R0_jikong2 and X12[i]>tschangdu]
    j=0
    for i in range(len(weizhi00)):
        if X12[weizhi00[i-j]]>tschangdu:
            del (X12[weizhi00[i-j]])
            del (R12[weizhi00[i-j]])
            del (alpha1[weizhi00[i-j]])
            del (h1[weizhi00[i-j]])
            j+=1
            # weizhi00=[weizhi00[i]-1 for i in range(len(weizhi00))]


    ######################################################################################################################



    # # 对厚度进行平稳处理
    # h1=smooth(h1)

    CPointX=[]
    CPointY=[]
    # 用B样条拟合外轮廓
    CPointX = X1[:]  # 整理B样条需要数据
    CPointY = R1[:]

    # 根据两边极孔条件删除不在范围内的点
    weizhi00=[i for i in range(len(CPointY)) if CPointY[i]<=R0_jikong1 and CPointX[i]<0]
    j=0
    for i in range(len(weizhi00)):
        del (CPointX[weizhi00[i-j]])
        del (CPointY[weizhi00[i-j]])
        j+=1
        # weizhi00=[weizhi00[i]-1 for i in range(len(weizhi00))]
    weizhi00=[i for i in range(len(CPointY)) if CPointY[i]<=R0_jikong2  and CPointX[i]>tschangdu]
    j=0
    for i in range(len(weizhi00)):
        if CPointX[weizhi00[i-j]]>tschangdu:
            del (CPointX[weizhi00[i-j]])
            del (CPointY[weizhi00[i-j]])
            j += 1
            # weizhi00=[weizhi00[i]-1 for i in range(len(weizhi00))]

    # 得到B样条的数据点
    knot = guifanhua(CPointX[:], CPointY[:])  # 生成B样条节点
    B3 = B3Spline(CPointX[:], CPointY[:], knot[:],R0_jikong2,Bsjindu)  # 得到B样条数据点

    # 2020.9.19改
    BsplineX=B3[0][:]
    BsplineY=B3[1][:]
    x1 = BsplineX[0] + (R0_jikong1 - BsplineY[1]) * ((BsplineX[0] - BsplineX[1]) / (BsplineY[0] - BsplineY[1]))
    BsplineX.insert(0,x1)
    BsplineY.insert(0,R0_jikong1)
    B3[0][:]=BsplineX
    B3[1][:]=BsplineY


    ######################################################################################################################


    # 内轮廓
    neisX=B3[0][:]
    neisR=B3[1][:]
    # plt.plot(neisX,neisR)
    # 端点
    # 左
    XddianZ.append(B3[0][0])
    RddianZ.append(B3[1][0])
    # 右
    XddianY.append(B3[0][-1])
    RddianY.append(B3[1][-1])


    ######################################################################################################################


    #外表面
    # 得到螺旋缠绕的外表面轮廓数据
    # 修正现在的X,R值，用B样条中的离该点X坐标最近的那个数据点代替原始X,R值
    cc = XiuZheng(X12, B3[0], B3[1])
    xiuzX = cc[0][:]
    xiuzR = cc[1][:]

    bb=quchong(xiuzX[:], xiuzR[:], h1[:],alpha1[:])#在计算出的内层B样条上找到当前反算外轮廓的数据点，去除重复值
    xiuzX =bb[0][:]
    xiuzR = bb[1][:]
    h1= bb[2][:]
    alpha1= bb[3][:]


    ######################################################################################################################


    # 得到环向层的厚度之后加入到螺旋层的筒身厚度计算中去，正角度的，负角度的后面计算自动加倍
    huanxsl=puceng.count('Hoop')
    tshouduhx=huanxsl*hxdt
    # 螺旋总厚度
    luoxuansl=puceng.count('HelixLoop')
    h2=[i*luoxuansl for i in h1]


    ######################################################################################################################


    # 后面需要处理的厚度,用来做赤道圆过渡处的处理
    XX11=[]
    RR11=[]
    # 找出筒身的位置
    tsX=[]
    tsR=[]
    XX =fengtouchazuo  ###左边赤道圆位置
    zuotongshen=0
    jdz = [abs(B3[0][j] - XX) for j in range(len(B3[0]))]  # 绝对值
    weizhi1=jdz.index(min(jdz))  # 获取各个x对应为修正位置

    #得到左部赤道圆位置坐标
    # tsX.append(B3[0][weizhi1])
    # tsR.append(B3[1][weizhi1])
    tsX.append(0)
    tsR.append(B3[1][weizhi1])

    # 提取筒身螺旋单层厚度/角度
    weizhi=[i for i in range(len(xiuzX)) if 0<xiuzX[i]<tschangdu]
    lxdt=h1[(weizhi[0]+weizhi[-1])/2]
    luoxuanjiao=abs(alpha1[(weizhi[0]+weizhi[-1])/2])


    # 左
    a=[xiuzX[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0]]
    b=[xiuzR[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0]]
    c=[h2[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0]]
    X=list(reversed(a))
    R=list(reversed(b))
    H1=list(reversed(c))
    waizuo=wailunkuozuo(X,R,H1,tsX)#得到极孔处修正的外轮廓，删掉了一部分厚度
    Xzuo=list(reversed(waizuo[0]))
    Rzuo=list(reversed(waizuo[1]))
    Hzuo=list(reversed(H1))
    j=0
    for i in range(len(Xzuo)):
        if Xzuo[i-j] <X[0]-fengtouchazuo-lxdt*luoxuansl*max(Rzuo)/R0_jikong1:#长度限定,对内外轮廓多做截断处理
            #外部
            del(Xzuo[i-j])
            del(Rzuo[i-j])
            del(Hzuo[i-j])
            del(a[i-j])
            del(b[i-j])
            j+=1
        if Xzuo[i-j]>tsX[0]:
            del (Xzuo[i - j])
            del (Rzuo[i - j])
            del (Hzuo[i - j])
            del (a[i - j])
            del (b[i - j])
            j += 1
    [XX11.append(a[i]) for i in range(len(a)) if (a[i]<=max(Xzuo)) and (a[i]>=min(Xzuo))]
    [RR11.append(b[i]) for i in range(len(a)) if (a[i]<=max(Xzuo)) and (a[i]>=min(Xzuo))]
    ############################################################
    # 得到右部赤道圆位置坐标
    youtongshen=zuotongshen+tschangdu
    jdz = [abs(B3[0][j] - XX) for j in range(len(B3[0]))]  # 绝对值
    weizhi2=jdz.index(min(jdz))  # 获取各个x对应为修正位置
    tsX.append(tsX[0] + tschangdu)
    tsR.append(tsR[0])#保证两侧赤道圆直径相等

    # 右
    X=[xiuzX[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    R=[xiuzR[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    Hyou=[h2[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    waiyou=wailunkuoyou(X,R,Hyou,tsX)#得到极孔处修正的外轮廓，删掉了一部分厚度
    Xyou=waiyou[0]
    Ryou=waiyou[1]
    j=0
    for i in range(len(Xyou)):
        if Xyou[i-j] >tsX[1]+fengtouchayou+lxdt*luoxuansl*max(Ryou)/R0_jikong2:#长度限定
           # 外部
            del (Xyou[i - j])
            del (Ryou[i - j])
            del (Hyou[i - j])
            del (X[i - j])
            del (R[i - j])
            j += 1
        if Xyou[i-j]<tsX[1]:
            del (Xyou[i - j])
            del (Ryou[i - j])
            del (Hyou[i - j])
            del (X[i - j])
            del (R[i - j])
            j += 1
    # 筒身数据
    [XX11.append(xiuzX[i]) for i in range(len(xiuzX)) if (xiuzX[i]>=tsX[0]) and (xiuzX[i]<=tsX[1])]
    [RR11.append(xiuzR[i]) for i in range(len(xiuzX)) if (xiuzX[i]>=tsX[0]) and (xiuzX[i]<=tsX[1])]
    # 右封头
    [XX11.append(X[i]) for i in range(len(X)) if (X[i]>=min(Xyou)) and (X[i]<=max(Xyou))]
    [RR11.append(R[i]) for i in range(len(X)) if (X[i]>=min(Xyou)) and (X[i]<=max(Xyou))]


    ######################################################################################################################


    # 删除掉极孔处有问题的数据点之后，重新组装新的外部轮廓数据
    #重新定义新的厚度数据，保护原始数据
    X11=[]
    R11=[]
    H11=[]

    #左封头数据
    [X11.append(i) for i in Xzuo]
    [R11.append(i) for i in Rzuo]
    [H11.append(i) for i in Hzuo]

    #筒身数据
    for i in range(len(xiuzX)):
        if xiuzX[i]<=tsX[1] and xiuzX[i]>=tsX[0]:
            X11.append(xiuzX[i])
            # R11.append(xiuzR[i])
            R11.append(tsR[0]+tshouduhx)
            # H11.append(h2[i])
            H11.append(lxdt*luoxuansl)

    #右封头数据
    [X11.append(i) for i in Xyou]
    [R11.append(i) for i in Ryou]
    [H11.append(i) for i in Hyou]


    ######################################################################################################################


    # 做赤道圆的过渡处理
    # 2020.9.22
    # 在赤道圆过渡位置寻找使其与封头相切的位置近似点
    leftTransitionX=tsX[0]
    leftTransitionR=tsR[0]+lxdt*luoxuansl*2+tshouduhx
    rightTransitionX=tsX[1]
    rightTransitionR=tsR[0]+lxdt*luoxuansl*2+tshouduhx
    tangentlineL=[]
    tangentlineR=[]
    wz=[]
    for i in range(len(X11)):
        if X11[i]<leftTransitionX:
            tangentlineL.append(math.atan(((X11[i] -leftTransitionX) / (R11[i] -leftTransitionR)))*180/math.pi)
        if X11[i] > rightTransitionX:
            tangentlineR.append(math.atan(((X11[i] -rightTransitionX)/(R11[i] - rightTransitionR)))*180/math.pi)
        if (X11[i] <= rightTransitionX) and (X11[i]>=leftTransitionX):
            wz.append(i)
    wz1=tangentlineL.index(max(tangentlineL))
    wz2=tangentlineR.index(min(tangentlineR))

    # 左
    for num in arange(wz1,min(wz)-1,1):
        t1=(tshouduhx)/(leftTransitionX-X11[wz1])*(X11[num]-X11[wz1])
        H11[num]=H11[num]+t1*0.5

    # 右
    for num in arange(max(wz)+1,max(wz)+1+wz2,1):
        t1=(tshouduhx)/(rightTransitionX-X11[max(wz)+1+wz2])*(X11[num]-X11[max(wz)+1+wz2])
        H11[num]=H11[num]+t1*0.5

    # 筒身部分环向层做加上环向层厚度的处理
    for num in arange(min(wz)-1,max(wz)+1,1):
        H11[num]=H11[num]+tshouduhx*0.5


    ######################################################################################################################


    # 对厚度光滑处理减缓厚度突变
    # Hh=smooth(H11)#光滑处理一下
    Hh=H11
    Xnei=[]
    Rnei=[]
    for ii in range(len(XX11)):#通过内部轮廓数据截断外部X两侧，后面做修正处理
        if XX11[ii]<=max(X11) and XX11[ii]>=min(X11):
        # if XX11[ii]<=tsX[1]+fengtouchayou and XX11[ii]>=tsX[0]-fengtouchazuo:#2020.9.29改
            Xnei.append(XX11[ii])
            Rnei.append(RR11[ii])
    try:
        wailunkuo = wailunkuo1(Xnei[:], Rnei[:], Hh[:])
    except:
        print("the wailunkuo command have problem")
        # 去除掉外轮廓x坐标在内轮廓横坐标控制范围外的点
    j=0
    for ii in range(len(wailunkuo[0][:])):
        if wailunkuo[0][ii-j]<min(Xnei[:]) or  wailunkuo[0][ii-j]>max(Xnei[:]):
            del(wailunkuo[0][ii-j])
            del(wailunkuo[1][ii-j])
            j += 1
    WX = wailunkuo[0][:]
    WR = wailunkuo[1][:]


    ######################################################################################################################


    # 外轮廓点的平稳处理
    # 左封头
    xzuo=[WX[i] for i in range(len(WX[:])) if WX[i]<=tsX[0]]
    rzuo=[WR[i] for i in range(len(WX[:])) if WX[i]<=tsX[0]]
    chuli=houdupingwen(xzuo[::-1],rzuo[::-1],pingwenxishu,tsX)
    WXzuo = chuli[0][::-1]
    WRzuo = chuli[1][::-1]
    # 筒身
    WXzhong=[WX[i] for i in range(len(WX[:])) if tsX[0]<WX[i]<tsX[1]]
    WRzhong=[WR[i] for i in range(len(WX[:])) if tsX[0]<WX[i]<tsX[1]]
    # 右封头
    xyou=[WX[i] for i in range(len(WX[:])) if WX[i]>=tsX[1]]
    ryou=[WR[i] for i in range(len(WX[:])) if WX[i]>=tsX[1]]
    chuli=houdupingwen(xyou[:],ryou[:],pingwenxishu,tsX)
    WX = WXzuo+WXzhong+chuli[0][:]
    WR = WRzuo+WRzhong+chuli[1][:]


    ######################################################################################################################


    # 这里的抛物线修正
    # 左侧
    zuoX11=min(neisX)-Hzuo[0]#lxdt*luoxuansl*max(WR)/R0_jikong1#左边极孔纤维堆积后的坐标
    # 最后两点的切线值
    Kk=(WX[0]-WX[1])/(WR[0]-WR[1])
    aa=(WX[0]-zuoX11-Kk*(WR[0]-R0_jikong1))/(((WR[0]**2-R0_jikong1**2))-2*WR[0]*(WR[0]-R0_jikong1))
    bb=Kk-2*aa*WR[0]
    cc=zuoX11-aa*R0_jikong1**2-bb*R0_jikong1
    XXjikongzuo=[]
    RRjikongzuo=[]
    for i in arange(R0_jikong1,WR[0],0.05):#靠近极孔这边的抛物线修正
        XXjikongzuo.append(aa*i**2+bb*i+cc)
        RRjikongzuo.append(i)
    WX=XXjikongzuo+WX
    WR=RRjikongzuo+WR

    # 右侧
    youX11=max(neisX)+Hyou[-1]#lxdt*luoxuansl*max(WR)/R0_jikong2#右边极孔纤维堆积后的坐标
    # 最后两点的切线值
    Kk=(WX[-2]-WX[-1])/(WR[-2]-WR[-1])
    aa=(WX[-1]-youX11-Kk*(WR[-1]-R0_jikong2))/(((WR[-1]**2-R0_jikong2**2))-2*WR[-1]*(WR[-1]-R0_jikong2))
    bb=Kk-2*aa*WR[-1]
    cc=youX11-aa*R0_jikong2**2-bb*R0_jikong2
    XXjikongyou=[]
    RRjikongyou=[]
    for i in arange(WR[-1],R0_jikong2,-0.05):
        XXjikongyou.append(aa*i**2+bb*i+cc)
        RRjikongyou.append(i)
    WX=WX+XXjikongyou
    WR=WR+RRjikongyou


    ######################################################################################################################


    # 对赤道圆处的赤道数据轮廓不明显采取的赤道圆加密数据处理，得到合适赤道圆数据
    # 2020.9.22
    CPointX = WX[:]
    CPointY =WR[:]
    # 筒身赤道处加密数据控制点，减少这里的变化突变
    # 找到赤道圆在控制点中的位置
    where1=[i for i in range(len(CPointX)-1) if CPointX[i]<tsX[0]<CPointX[i+1]]#左赤道圆
    where2=[i for i in range(len(CPointX)-1) if CPointX[i]<tsX[1]<CPointX[i+1]]#右赤道圆
    CPointX.insert(where1[0]+1,tsX[0])#扩充赤道圆值
    CPointY.insert(where1[0]+1,CPointY[where1[0]+1])#扩充赤道圆值
    CPointX.insert(where2[0]+2,tsX[1])
    CPointY.insert(where2[0]+2,CPointY[where2[0]+1])#扩充赤道圆值
    where1=[i for i in range(len(CPointX)) if CPointX[i]==tsX[0]]#左赤道圆
    where2=[i for i in range(len(CPointX)) if CPointX[i]==tsX[1]]#右赤道圆


    ######################################################################################################################


    # 赤道圆过度处加密
    p=5
    qujian=[where1[0],where2[0]]
    for i in qujian:
        if i==qujian[0]:
            for j in arange(1,p,1):
                [CPointX.insert((i+1),CPointX[(i)]+j*(CPointX[(i+1)]-CPointX[(i)])/p)]
                [CPointY.insert((i+1),CPointY[(i+1)])]
        if i == qujian[1]:
            for j in arange(1, p, 1):
                [CPointX.insert((i +p-1), CPointX[(i +p-1)]-j * (CPointX[(i + p-1)] - CPointX[(i+p-2)]) / p)]
                [CPointY.insert((i +p-1), CPointY[(i - 2+p)])]


    ######################################################################################################################


    # 生成外轮廓B样条节点
    knot = guifanhua(CPointX[:], CPointY[:])  # 生成B样条节点
    B3 = B3Spline(CPointX[:], CPointY[:], knot[:],R0_jikong2,Bsjindu)  # 得到B样条数据点
    # 2020.9.19改
    BsplineX=B3[0][:]
    BsplineY=B3[1][:]
    x1 = BsplineX[0] + (R0_jikong1 - BsplineY[1]) * ((BsplineX[0] - BsplineX[1]) / (BsplineY[0] - BsplineY[1]))
    BsplineX.insert(0,x1)
    BsplineY.insert(0,R0_jikong1)
    B3[0][:]=BsplineX
    B3[1][:]=BsplineY


    ######################################################################################################################


    # 得到外轮廓绘制图形数据
    waisX=B3[0][:]
    waisR=B3[1][:]
    # plt.plot(waisX,waisR)
    # 端点
    # 左
    XddianZ.append(B3[0][0])
    RddianZ.append(B3[1][0])
    # 右
    XddianY.append(B3[0][-1])
    RddianY.append(B3[1][-1])
    # 绘图

    # 从B样条数据里面取一部分，防止导入abaqus卡死

    neisXR=Qujiange(neisX,neisR,len(X12))
    neisX=neisXR[0][:]
    neisR=neisXR[1][:]

    waisXR=Qujiange(waisX,waisR,len(X12))
    waisX=waisXR[0][:]
    waisR=waisXR[1][:]

    # plt.plot(neisX,neisR,'*y')
    # plt.plot(waisX,waisR,'*b')
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()

    # plt.show()


    ######################################################################################################################


    # 二、下面是abaqus加墨分析模块
        
        
    ######################################################################################################################
    # abaqus建模
    #引入abaqus模块建模
    # #导入abaqus模块
    # from abaqus import *
    # from abaqusConstants import *
    # from caeModules import *


    Mdb()#创建新的模型数据库
    modelName = 'composite'
    m = mdb.Model(name=modelName)


    #在abaqus中绘制草图
    session.viewports['Viewport: 1'].view.fitView()
    s1 = m.ConstrainedSketch(name='composite-'+str('PressureVessel'), sheetSize=int(waisX[-1]))
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    constructionline=s1.ConstructionLine(point1=(waisX[0], 0), point2=(waisX[-1], 0))  # 根据外层X坐标自动调整轴线长度
    s1.assignCenterline(line=constructionline)
    s1.FixedConstraint(entity=constructionline)
    Vline=s1.ConstructionLine(point1=(0.0, -20), point2=(0.0, max(waisR)))
    s1.FixedConstraint(entity=Vline)


    ######################################################################################################################




    # j=0
    # for i in range(1,len(neisX)-2):
    #     if waisX[i+1-j]<tsX[0] and waisX[i-j]>min(Xnei[:]):
    #         if (waisX[i+1-j]-waisX[i-j])**2+(waisR[i+1-j]-waisR[i-j])**2<1:#距离太近的清除掉,防止后面abaqus的B样条畸形
    #             del(waisX[i-j])
    #             del(waisR[i-j])
    #             j += 1
    #     if waisX[i+1-j]>tsX[1] and waisX[i+1-j]<max(Xnei[:]):
    #         if (waisX[i+1-j]-waisX[i-j])**2+(waisR[i+1-j]-waisR[i-j])**2<1:
    #             del(waisX[i+1-j])
    #             del(waisR[i+1-j])
    #             j += 1


    ######################################################################################################################


    Spoints = zip(neisX[:], neisR[:])#打包坐标点
    curve1 = s1.Spline(points=Spoints)
    Spoints = zip(waisX[:], waisR[:])#打包坐标点
    curve2 = s1.Spline(points=Spoints)
    line_cy5 = s1.Line(point1=(XddianZ[0], RddianZ[0]),
                       point2=(XddianZ[1], RddianZ[1]))  # 端点的存储没有区分纯螺旋与非纯螺旋
    line_cy6 = s1.Line(point1=(XddianY[0], RddianY[0]),
                       point2=(XddianY[1], RddianY[1]))

    # 生成三维实体模型
    session.viewports['Viewport: 1'].view.fitView()
    s1.setPrimaryObject(option=STANDALONE)
    s1.sketchOptions.setValues(constructionGeometry=ON)
    p = mdb.models['composite'].Part(name='Part-'+str('PressureVessel'), dimensionality=THREE_D,type=DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch=s1, angle=zhuanjiao, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()


    ######################################################################################################################


    # 导入装配模块
    a = mdb.models['composite'].rootAssembly
    #session.viewports['Viewport: 1'].setValues(displayedObject=a)#显示，#转换到装配体模块视图
    #session.viewports['Viewport: 1'].assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['composite'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)


    ######################################################################################################################


    # 调节显示精度
    p = mdb.models['composite'].parts['Part-'+str('PressureVessel')]
    p.setValues(geometryRefinement=EXTRA_FINE)
    a.Instance(name='Part-'+str('PressureVessel')+str(1), part=p, dependent=ON)


    ######################################################################################################################

    # 定义叠层上表面
    f,e=p.faces,p.edges
    lenn=len(waisX)/2
    lenn1=len(neisX)/2
    side1Faces22 = f.findAt(((waisX[lenn],waisR[lenn]*cos(0.5*zhuanjiao*todegree),waisR[lenn]*sin(0.5*zhuanjiao*todegree)), ))
    side1Faces22=side1Faces22+f.findAt(((waisX[lenn/2],waisR[lenn/2]*cos(0.5*zhuanjiao*todegree),waisR[lenn/2]*sin(0.5*zhuanjiao*todegree)), ))
    side1Faces22=side1Faces22+f.findAt(((waisX[4*lenn/3],waisR[4*lenn/3]*cos(0.5*zhuanjiao*todegree),waisR[4*lenn/3]*sin(0.5*zhuanjiao*todegree)), ))
    # for ii in range(1:len(waisX)):
    p.Surface(side1Faces=side1Faces22, name='Surf-top')#后面赋予属性可以用


    ######################################################################################################################


    # 切片疏密大小控制
    zuofengtou=[]
    youfengtou=[]
    juli=1
    yuanshishuju=zip(X12,R12,alpha1,h2)#得到赋予角度的原始cadwind数据
    pingjudistance=0.5#得到这里的原始数据点平均距离
    shanchushu=[]
    for i in range(1,len(yuanshishuju)):
        distan=((yuanshishuju[i][0]-yuanshishuju[i-1][0])**2+(yuanshishuju[i][1]-yuanshishuju[i-1][1])**2)**0.5
        if distan<pingjudistance:#片区切片大小控制
            shanchushu.append(yuanshishuju[i-1])
    # 距离较近两点取舍，删除距离较近的点
    for i in shanchushu[::-1]:
        yuanshishuju.remove(i)

    # 根据赤道圆数据坐标，切出左右封头数据
    for i in range(len(yuanshishuju)):
        if yuanshishuju[i][0]<=zuotongshen :
            # zuofengtou=yuanshishuju[0:i+1]
            zuofengtou.append(yuanshishuju[i])
            # break
    for i in range(len(yuanshishuju)):
    # for i in range(len(X12)):
        if yuanshishuju[i][0]>=youtongshen:
            # youfengtou=yuanshishuju[i:-1]
            youfengtou.append(yuanshishuju[i])
            # break


    ######################################################################################################################


    # 创建切分直线，为切分封头实体模型做准备
    def pressureqieFen(fengtou):
        top_xvalue,top_yvalue,bot_xvalue,bot_yvalue=[],[],[],[]
        for i in range(len(fengtou)-1):
            k1=(fengtou[i+1][1]-fengtou[i][1])/(fengtou[i+1][0]-fengtou[i][0])
            if k1==0:
                top_xvalue.append(fengtou[i][0])
                top_yvalue.append((fengtou[i][1]+(huanxsl+luoxuansl)*fengtou[i][3]))
                bot_xvalue.append(fengtou[i][0])
                bot_yvalue.append((fengtou[i][1]-(huanxsl+luoxuansl)*fengtou[i][3]))
            else:
                k2=-1/k1
                if fengtou[i][0]<=0 and i<(len(fengtou)-2):
                    theta1=arctan(k2)#弧度值
                    # print(theta1)
                    top_xvalue.append((fengtou[i+1][0]+(huanxsl+luoxuansl)*fengtou[i][3]*cos(theta1)))
                    top_yvalue.append((fengtou[i+1][1]+(huanxsl+luoxuansl)*fengtou[i][3]*sin(theta1)))
                    bot_xvalue.append((fengtou[i+1][0]-(huanxsl+luoxuansl)*fengtou[i][3]*cos(theta1)))
                    bot_yvalue.append((fengtou[i+1][1]-(huanxsl+luoxuansl)*fengtou[i][3]*sin(theta1)))
                if fengtou[i][0]<=0 and i==(len(fengtou)-2):
                    theta1 = arctan(k2)  # 弧度值
                    # print(theta1)
                    top_xvalue.append(fengtou[i + 1][0])
                    top_yvalue.append(fengtou[i + 1][1]+ (huanxsl + luoxuansl) * fengtou[i][3])
                    bot_xvalue.append(fengtou[i + 1][0])
                    bot_yvalue.append(fengtou[i + 1][1] - (huanxsl + luoxuansl) * fengtou[i][3])
                if fengtou[i][0] >= youtongshen and i==0:
                    theta1 = arctan(k2)  # 弧度值
                    top_xvalue.append(fengtou[i][0])
                    top_yvalue.append(fengtou[i][1] +(huanxsl+luoxuansl)*fengtou[i][3])
                    bot_xvalue.append(fengtou[i][0])
                    bot_yvalue.append(fengtou[i][1] -(huanxsl+luoxuansl)*fengtou[i][3])
                if fengtou[i][0] >= youtongshen and i >=1:
                    theta1 = arctan(k2)  # 弧度值
                    top_xvalue.append((fengtou[i][0] +(huanxsl+luoxuansl)*fengtou[i][3] * cos(theta1)))
                    top_yvalue.append((fengtou[i][1] +(huanxsl+luoxuansl)*fengtou[i][3] * sin(theta1)))
                    bot_xvalue.append((fengtou[i][0] -(huanxsl+luoxuansl)*fengtou[i][3] * cos(theta1)))
                    bot_yvalue.append((fengtou[i][1] -(huanxsl+luoxuansl)*fengtou[i][3] * sin(theta1)))
        qiefenwailunkuo=zip(top_xvalue,top_yvalue)
        qiefenneilunkuo=zip(bot_xvalue,bot_yvalue)
        return qiefenwailunkuo,qiefenneilunkuo

    # 用厚度来分割
    zqiefenwailunkuo,zqiefenneilunkuo=pressureqieFen(zuofengtou)
    yqiefenwailunkuo,yqiefenneilunkuo=pressureqieFen(youfengtou)


    ######################################################################################################################


    # 进行切分侧面处理处理
    f, e, dpart = p.faces, p.edges, p.datums

    # 建立分割草图定位面
    pickFace=f.findAt(coordinates=((0.5*(tsX[0]+tsX[1]),0.5*(max(neisR)+max(waisR)), 0.0)))

    s = mdb.models['composite'].ConstrainedSketch(name='PartitionSketch',
        sheetSize=(max(waisX)-min(waisX))*3, gridSpacing=round((max(waisX)-min(waisX))/10))
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    # 修建分割草图
    try:
        for k in range(0,len(zqiefenwailunkuo)):#建立分割草图线段
            line1=s.Line(point1=zqiefenwailunkuo[k],point2=zqiefenneilunkuo[k])
            s.autoTrimCurve(curve1=line1, point1=zqiefenwailunkuo[k])#修建直线
            line11=g.findAt(zqiefenneilunkuo[k])#寻找修建后的直线
            s.autoTrimCurve(curve1=line11, point1=zqiefenneilunkuo[k])
        for k in range(0,len(yqiefenwailunkuo)):#建立分割草图线段
            line2=s.Line(point1=yqiefenwailunkuo[k],point2=yqiefenneilunkuo[k])
            s.autoTrimCurve(curve1=line2, point1=yqiefenwailunkuo[k])
            line22 = g.findAt(yqiefenneilunkuo[k])  # 寻找修建后的直线
            s.autoTrimCurve(curve1=line22, point1=yqiefenneilunkuo[k])
    except:
        print("Partition Error")
    # 定义一条边用于草图定位显示
    Axisline=p.DatumAxisByTwoPoint(point1=(waisX[0],0,0),  point2=(waisX[-1],0,0))#创建轴线
    p.PartitionFaceBySketch(sketchUpEdge=dpart[Axisline.id],faces=pickFace,sketchOrientation=TOP, sketch=s)


    ######################################################################################################################


    # 切分模型单元
    qiefenbian=[]
    quxianbian=[]

    f1, e1, d2 = p.faces, p.edges, p.datums
    for edge in e1:
        try:
            if abs(edge.pointOn[0][2])<1e-6 and edge.pointOn[0][1]>0 :#找到z=0和半径大于0的面，z有偏差，
                qiefenbian.append(edge)
                edge.getCurvature(parameter=0.5)#尝试有曲率的边，0.5为尝试数，只要不是0
                quxianbian.append(edge)
        except:
            pass
    quchu=[qiefenbian.remove(x) for x in quxianbian]#去除曲线边

    refpoint=[]
    for i in range(len(qiefenbian)):
        dian=qiefenbian[i].pointOn
        refpoint.append(dian[0])#获得切分边直线边上的第一号数据点

    # 切分单元
    for point in refpoint:
        try:#尝试分割实体模型
            c=p.cells               
            pickedEdges =(e1.findAt(coordinates=point),)   
            p.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(neisX[-1], neisR[-1]*cos(0.5*zhuanjiao*todegree), neisR[-1]*sin(0.5*zhuanjiao*todegree))), cells=c, edges=pickedEdges)
        except:#尝试切分模型，遇到错误信息
            pass


    ######################################################################################################################


    # 为了在筒身中部添加约束，保证两边极孔的受力形式一样，需将边界约束添加到筒身中部
    point=(youtongshen/2, 0.0, 0.0)
    Planepoint=p.DatumPointByCoordinate(coords=point)#创建过平面的点
    Axisline=p.DatumAxisByTwoPoint(point1=(waisX[0],0,0),  point2=(waisX[-1],0,0))#创建轴线
    d_part = p.datums
    PlaneByPoint=p.DatumPlaneByPointNormal(point=d_part[Planepoint.id], normal=d_part[Axisline.id])
    # 分割中间cell
    pickedCells = c.findAt((youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
    , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180), ))#得到筒身cell
    p.PartitionCellByDatumPlane(datumPlane=d_part[PlaneByPoint.id], cells=pickedCells)


    ######################################################################################################################


    # 得到各个小片区的中点坐标
    # 初始化相关的曲面数据，后面需要
    Masterface=[]
    # loadsurfacecenterpoint=[]
    f_all=p.faces#得到part所有的面
    face_coordinate=[]
    face_coordinate1=[]
    face_coordinate2=[]
    face_coordinate3=[]

    for face in f_all:#寻找每一个面的中心点
        centercoordinate=face.getCentroid()#将曲面一个一个的输入判断
        # print (centercoordinate)
        if centercoordinate[0][2] ==0.0 and centercoordinate[0][1]>0.0:#根据面的中心点，得到z为0的面
            Masterface.append(centercoordinate)
            if centercoordinate[0][0]<tsX[0]:
                face_coordinate1.append(centercoordinate[0])
            if centercoordinate[0][0]>=tsX[0] and centercoordinate[0][0]<=tsX[1]:
                face_coordinate2.append(centercoordinate[0])
            if centercoordinate[0][0]>tsX[1]:
                face_coordinate3.append(centercoordinate[0])
    # face_coordinate1.sort(key=(lambda x: x[0]))#这里的角度排序有问题，因为X不是单调递增的，改动如下
    face_coordinate1.sort(key=(lambda x: x[1]))#这里的角度排序有问题，因为X不是单调递增的，封头按半径排序
    face_coordinate2.sort(key=(lambda x: x[0]))#这里的角度排序有问题，因为X不是单调递增的，筒身按X排序
    face_coordinate3.sort(key=(lambda x: -x[1]))#这里的角度排序有问题，因为X不是单调递增的，封头按半径排序
    # print(len(face_coordinate1))
    # print(face_coordinate1)
    # print(len(face_coordinate2))
    # print(face_coordinate2)
    # print(len(face_coordinate3))
    # print(face_coordinate3)
    for i in range(len(face_coordinate1)+len(face_coordinate2)+len(face_coordinate3)):
        if i <len(face_coordinate1):
            face_coordinate.append(face_coordinate1[i])
        if i < (len(face_coordinate1)+len(face_coordinate2)) and i >=len(face_coordinate1):
            face_coordinate.append(face_coordinate2[i-((len(face_coordinate1)+len(face_coordinate2)))])
        if i >= (len(face_coordinate1)+len(face_coordinate2)) and i<(len(face_coordinate1)+len(face_coordinate2)+len(face_coordinate3)):
            face_coordinate.append(face_coordinate3[i-((len(face_coordinate1)+len(face_coordinate2)))])
    # print(face_coordinate)


    ######################################################################################################################


    #创建复合材料工程常数
    mat_comp=m.Material(name='Material-composite')
    mat_comp.Elastic(type=ENGINEERING_CONSTANTS, 
        table=(engineeringdata, ))
    m.HomogeneousSolidSection(name='Section-composite', material='Material-composite', thickness=None)

    list_angle=[]
    for i in range(len(zuofengtou)):
        list_angle.append(zuofengtou[i][2])
    list_angle.append(0)
    for i in range(len(youfengtou)):
        list_angle.append(youfengtou[i][2])

    count=0
    tongshendian=[]#分别建立筒身和封头的set集
    for k in face_coordinate:
        c=p.cells      
        count=count+1
        pickcells=c.findAt(((k[0], k[1], k[2]),))
        p.Set(name='Set'+'-'+str(count),cells=pickcells)
        if count==len(zuofengtou)+1 or  count==len(zuofengtou)+2:
            tongshendian.append((k[0], k[1], k[2]),)
            

    ######################################################################################################################


    # 封头赋予铺层角度属性
    for k in range(1,count+1):
        if k!=len(face_coordinate1)+1:#只要不等于筒身集合就赋予封头属性
            layupOrientation = None
            region1 = p.sets['Set-'+str(k)]
            normalAxisRegion = p.surfaces['Surf-top']
            compositeLayup = p.CompositeLayup(
                name='CompositeLayup-'+str(k), description='', elementType=SOLID, 
                symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON,
                                   poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT,
                                   useDensity=OFF)
            for i in range(1,luoxuansl+1):
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*i-1), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0,
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1],
                    additionalRotationType=ROTATION_NONE, additionalRotationField='',
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*i), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0,
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1],
                    additionalRotationType=ROTATION_NONE, additionalRotationField='',
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None,
                additionalRotationType=ROTATION_NONE, angle=0.0,
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3,
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion,
                normalAxisDirection=AXIS_3, flipNormalDirection=False,
                primaryAxisDefinition=VECTOR, primaryAxisVector=(1, 0.0, 0.0),
               primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)
        else:
            pass


    ######################################################################################################################


    ###筒身赋予铺层角度属性
    layupOrientation = None
    for i in range(1,3):
        region1 = p.sets['Set-'+str(len(face_coordinate1)+i)]
        normalAxisRegion = p.surfaces['Surf-top']
        compositeLayup = p.CompositeLayup(
            name='CompositeLayup-'+str(len(face_coordinate1)+i), description='', elementType=SOLID,
            symmetric=False)
        compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
            poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
            useDensity=OFF)
        ply_count=0
        for x in puceng:
            if x=='HelixLoop':
                ply_count=ply_count+1
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*ply_count-1), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=lxdt, 
                    orientationType=SPECIFY_ORIENT, orientationValue=luoxuanjiao, 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-' + str(2*ply_count), region=region1,
                                            material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=lxdt,
                                            orientationType=SPECIFY_ORIENT, orientationValue=-luoxuanjiao,
                                            additionalRotationType=ROTATION_NONE, additionalRotationField='',
                                            axis=AXIS_3, angle=0, numIntPoints=3)
            elif x=='Hoop':
                ply_count=ply_count+1
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*ply_count-1), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=hxdt, 
                    orientationType=SPECIFY_ORIENT, orientationValue=90,
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0, numIntPoints=3)

        compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
            additionalRotationType=ROTATION_NONE, angle=0.0, 
            additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
            normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
            normalAxisDirection=AXIS_3, flipNormalDirection=False, 
            primaryAxisDefinition=VECTOR, primaryAxisVector=(1, 0.0, 0.0), 
           primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 


    ######################################################################################################################


    # 创建step
    mdb.models['composite'].StaticStep(name='Step-1', previous='Initial',
        maxNumInc=300, initialInc=0.05, maxInc=0.1)
    mdb.models['composite'].steps['Step-1'].setValues(nlgeom=ON)

    ############################################################
    # 创建场变量输出
    for i in range(2,len(face_coordinate)+2):
        mdb.models['composite'].FieldOutputRequest(name='F-Output-'+str(i),
            createStepName='Step-1', variables=('S', 'E', 'LE'), layupNames=(
            'Part-PressureVessel1.CompositeLayup-'+str(i-1), ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True,
            outputAtPlyBottom=False, rebar=EXCLUDE)


    ######################################################################################################################


    # 通过旋转侧面的中点坐标创建柱坐标
    # DatumcsysCylinder=a.DatumCsysByThreePoints(\
    #     point1=(face_coordinate2[0][0],face_coordinate2[0][1]*math.cos(zhuanjiao/2*math.pi/180),\
    #                                 face_coordinate2[0][1]*math.sin(zhuanjiao/2*math.pi/180)), name='Datum csys-Cylinder',
    #     coordSysType=CYLINDRICAL, origin=(face_coordinate2[0][0], 0.0, 0.0),
    #     point2=face_coordinate2[0])

    DatumcsysCylinder=a.DatumCsysByThreePoints(\
        point1=(youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
        , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180)), name='Datum csys-Cylinder1',
        coordSysType=CYLINDRICAL, origin=(youtongshen/2, 0.0, 0.0),
        point2=(youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
        , 0))


    ######################################################################################################################


    # 创建参考点
    rp=a.ReferencePoint(point=point)
    ReferencePoints=a.referencePoints#assembly层次的参考点存储仓库
    region1=a.Set(referencePoints=[ReferencePoints[rp.id],], name='m_Set-1')
    f1 = a.instances['Part-PressureVessel1'].faces
    faceszhong = f1.findAt(((youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
    , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180)), ))


    ######################################################################################################################


    # 找到主面的中点坐标
    Slaveface=[]
    for i in range(len(Masterface)):
        abc1=Masterface[i][0][0]
        abc2=Masterface[i][0][1]*math.cos(math.pi*zhuanjiao/180)
        abc3=Masterface[i][0][1]*math.sin(math.pi*zhuanjiao/180)
        Slaveface.append(((abc1,abc2,abc3), ))
    # 创建主从面
    s1=a.instances['Part-PressureVessel1'].faces
    side1Faces1=s1.findAt(Masterface[0])
    for i in range(1,len(Masterface)):
        faces1=s1.findAt(Masterface[i])
        side1Faces1=side1Faces1+faces1#遍历读取面
    regionMaterface=a.Surface(side1Faces=side1Faces1, name='m_Surf-1')
    side1Faces2=s1.findAt(Slaveface[0])
    for i in range(1,len(Slaveface)):
        faces2 = s1.findAt(Slaveface[i])
        side1Faces2+=faces2#遍历读取面
    regionSlaveface=a.Surface(side2Faces=side1Faces2, name='s_Surf-1')
    # 创建参考带点，循环对称用
    r1 = a.ReferencePoint(point=(0.0, 0.0, 0.0))
    r1points=a.referencePoints
    regionreferencepoint1=regionToolset.Region(referencePoints=[r1points[r1.id]])#该种做法参见《abaquspython二次开发攻略》203页
    # region3=regionToolset.Region(referencePoints=refPoints1)
    r1 = a.ReferencePoint(point=(100.0, 0.0, 0.0))
    rpoints=a.referencePoints
    # refPoints1=(r1[2], )
    regionreferencepoint2=regionToolset.Region(referencePoints=[r1points[r1.id]])#该种做法参见《abaquspython二次开发攻略》203页
    # region4=regionToolset.Region(referencePoints=refPoints1)
    # 建立循环对称相互作用
    mdb.models['composite'].CyclicSymmetry(name='Int-1', createStepName='Initial',
        master=regionMaterface, slave=regionSlaveface, axisPoint1=regionreferencepoint1, axisPoint2=regionreferencepoint2,
        positionToleranceMethod=COMPUTED_TOLERANCE, positionTolerance=0.0,
        adjustTie=True, repetitiveSectors=360/int(zhuanjiao),
        extractedNodalDiameter=ALL_NODAL_DIAMETER, excitationNodalDiameter=0)
    #: The interaction "Int-1" has been created.


    ######################################################################################################################


    # 由前面的内轮廓数据点坐标，得到内轮廓曲面旋转后的每个小面的中点坐标
    # 将建立内表面的
    # 2020.9.26改，给每个cell找到面，筛选出内外表面，进而得到内表面
    # 1.遍历每一个cell
        ##1.1得到cell的每1个面
        ##1.2根据曲面的曲率判定上下弧面,内外弧面含有两个曲率半径，其他要么是平面要么其他曲面
    #2.将得到的内外表面上数据点按径向距离区分内外表面
    # 建立内表面的几何surfaces
    neiwaiSurfaces=[]
    neisurfaces=[]
    for iijj  in range((len(face_coordinate))):#cell的数量是和侧面中点数是一样的
        cellface=p.cells[iijj].getFaces()#得到当前cell的6个表面
        facex = []
        neiwai=[]
        for ih in cellface:
            facepoint0 = p.faces[ih].pointOn[0]  # 得到当前cell当前编号face的随机面内坐标,不能用中点，因为曲面的中点不一定在面上
            try:
                p.faces[ih].getCurvature(0, 0)["curvature1"] != 0 or p.faces[ih].getCurvature(0, 0)["curvature2"] != 0
                qulv1 = p.faces[ih].getCurvature(0, 0)["curvature1"]
                # print(qulv1)
                qulv2 = p.faces[ih].getCurvature(0, 0)["curvature2"]
                # print(qulv2)
                if qulv1  != 0 and qulv2 != 0:
                    neiwai.append((facepoint0,ih))
                    facex.append(p.faces[ih].pointOn[0][0])
                    # print(neiwai)
            except:
                pass
        if len(facex)==2:
            if (facex[0] + facex[1]) / 2 >= 0 and (facex[0] + facex[1]) / 2 <= youtongshen:  # 先判断是不是筒身，筒身的中点连判断，因为筒身这里的判断有嗲小问题不用这个的话
                print(facex)
                neiwai1=[neiwai[i][0] for i in range(2)]
                neiwai1.sort(key=(lambda x: (x[1] ** 2 + x[2] ** 2)))  # 按旋转半径排序
                print(neiwai)
                print(neiwai1)
                # neisurfaces.append((neiwai1[0][0],))
                # neiwaiSurfaces.append(neiwai)  # 存储内外表面
    StackDirectionSurface = neiwai1[-1]
    neibiaomian = p.faces[neiwai[0][1]]
    neibiaomian1 = neibiaomian.getFacesByFaceAngle(10)


    #通过曲率识别
    # neiwaiSurfaces=[]
    # neisurfaces=[]
    # for iijj  in range((len(face_coordinate))):#cell的数量是和侧面中点数是一样的
    #     cellface=p.cells[iijj].getFaces()#得到当前cell的6个表面
    #     print(p.cells[iijj])
    #     print(iijj)
    #     facepoint0=[]#每个面1点一组存储
    #     facex=[]
    #     neiwai=[]
    #     for ih in cellface:
    #         facepoint0 = p.faces[ih].pointOn[0]  # 得到当前cell当前编号face的随机面内坐标,不能用中点，因为曲面的中点不一定在面上
    #         try:
    #             p.faces[ih].getCurvature(0.5, 0.5)["curvature1"] != 0 or p.faces[ih].getCurvature(0.5, 0.5)["curvature2"] != 0
    #             qulv1 = p.faces[ih].getCurvature(0.5, 0.5)["curvature1"]
    #             # print(qulv1)
    #             qulv2 = p.faces[ih].getCurvature(0.5, 0.5)["curvature2"]
    #             # print(qulv2)
    #             if qulv1  != 0 and qulv2 != 0:
    #                 neiwai.append(facepoint0)
    #                 facex.append(p.faces[ih].pointOn[0][0])
    #                 # print(neiwai)
    #         except:
    #             pass
    #             # print("平面")
    #     if len(facex)==2:
    #             # 左
    #         if (facex[0] + facex[1]) / 2 < zuotongshen:  # 先判断是不是左封头
    #             neiwai.sort(key=(lambda x: x[0] ** 2 + x[1] ** 2 + x[2] ** 2))  # 按旋转半径排序
    #             neisurfaces.append((neiwai[0],))
    #             print(neiwai)
    #             neiwaiSurfaces.append(neiwai)  # 存储内外表面
    #         # 筒身
    #         if (facex[0] + facex[1]) / 2 >= 0 and (facex[0] + facex[1]) / 2 <= youtongshen:  # 先判断是不是筒身，筒身的中点连判断，因为筒身这里的判断有嗲小问题不用这个的话
    #             neiwai.sort(key=(lambda x: (x[1] ** 2 + x[2] ** 2)))  # 按旋转半径排序
    #             neisurfaces.append((neiwai[0],))
    #             neiwaiSurfaces.append(neiwai)  # 存储内外表面
    #             StackDirectionSurface = neiwai[-1]
    #         # 右封头
    #         if (facex[0] + facex[1]) / 2 > youtongshen:  # 先判断是不是右封头
    #             neiwai.sort(key=(lambda x: (x[0] - youtongshen) ** 2 + x[1] ** 2 + x[2] ** 2))  # 按旋转半径排序
    #             neisurfaces.append((neiwai[0],))
    #             neiwaiSurfaces.append(neiwai)  # 存储内外表面
    #
    #     facex=[]
    #     neiwai=[]
    #     for ih in cellface:
    #         facepoint0 = p.faces[ih].pointOn[0]  # 得到当前cell当前编号face的随机面内坐标,不能用中点，因为曲面的中点不一定在面上
    #         try:
    #             p.faces[ih].getCurvature(0.5, 0.5)["curvature1"] == 0 and p.faces[ih].getCurvature(0.5, 0.5)["curvature2"] != 0
    #             qulv1 = p.faces[ih].getCurvature(0.5, 0.5)["curvature1"]
    #             # print(qulv1)
    #             qulv2 = p.faces[ih].getCurvature(0.5, 0.5)["curvature2"]
    #             # print(qulv2)
    #             # if qulv1 == 0:
    #             neiwai.append(facepoint0)
    #             facex.append(p.faces[ih].pointOn[0][0])
    #                 # print(neiwai)
    #         except:
    #             pass
    #             # print("平面")
    #     if len(facex) == 2:
    #         # 筒身
    #         if (facex[0] + facex[1]) / 2 >= 0 and (facex[0] + facex[1]) / 2 <= youtongshen:  # 先判断是不是筒身，筒身的中点连判断，因为筒身这里的判断有嗲小问题不用这个的话
    #             neiwai.sort(key=(lambda x: (x[1] ** 2 + x[2] ** 2)))  # 按旋转半径排序
    #             neisurfaces.append((neiwai[0],))
    #             neiwaiSurfaces.append(neiwai)  # 存储内外表面
    #             StackDirectionSurface = neiwai[-1]
    # side1Faces3= s1.findAt(neisurfaces[0])#通过中面寻找到曲面
    # for i in range(1,len(neisurfaces)):
    #     faces3=s1.findAt(neisurfaces[i])
    #     side1Faces3=side1Faces3+faces3#遍历读取面

    # session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry= COORDINATE)

    # region5=a.Surface(side1Faces=side1Faces3, name='Surf-bottom')

    p.Surface(side1Faces=neibiaomian1, name='Surf-bottom')
    region5=a.instances['Part-PressureVessel1'].surfaces['Surf-bottom']
    # mdb.models['composite'].loads['Load-1'].setValues(region=region5)
    mdb.models['composite'].Pressure(name='Load-1', createStepName='Step-1',
        region=region5, distributionType=UNIFORM, field='', magnitude=pressure,
        amplitude=UNSET)


    ######################################################################################################################


    # 建立边界条件1,对旋转周期对称约束进行边界约束
    regionBC = regionToolset.Region(faces=side1Faces1)
    datum = a.datums[DatumcsysCylinder.id]
    mdb.models['composite'].YsymmBC(name='BC-1', createStepName='Initial',
        region=regionBC, localCsys=datum)

    # 边界条件2，对极孔这里的的轴向进行限定
    regiontshen= regionToolset.Region(faces=faceszhong)
    mdb.models['composite'].DisplacementBC(name='BC-2', createStepName='Initial',
        region=regiontshen, u1=UNSET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET,
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)


    ######################################################################################################################


    # 种子点分布控制，厚度方向一个网格
    p.deleteMesh()
    # 思路：
    # 1.由之前旋转侧面的中点数量坐标得到cell数量
    # 2.得到每一个cell的边
    # 3.找到cell的各边单独控制网格数量
    # 建立cell的边数据
    cellzuoyouzong=[]#所有单元左右边
    cellshangxiazong=[]#所有单元左右边
    cellcurveszong=[]#所有单元左右边
    for iijj  in range((len(face_coordinate))):#cell的数量是和侧面中点数是一样的
        celledges=p.cells[iijj].getEdges()#得到当前cell的12个边
        edgepoint0=[]
        for ih in celledges:
            edgepoint0.append((p.edges[ih].pointOn[0],ih))#得到当前cell当前编号edge的线上随机坐标
    # 2020.10.22
    #     根据Z坐标排序
        edgepoint0.sort(key=(lambda x: x[0][2]))  #一个cell有12个边
        Z00=edgepoint0[0:4]
        ZFei0=edgepoint0[8::]
        ZHubian=edgepoint0[4:8]
        Z0=Z00+ZFei0
        Xx0=sum([Z0[i][0][0] for i in range(len(Z0))])/len(Z0)
        if tsX[0] < Xx0 < tsX[1]:  # 筒身
            Z0.sort(key=(lambda x: (x[0][1] ** 2 + x[0][2] ** 2)))  # 旋转半径排序
            # print(Z0)
            partedges = p.edges
            for i in range(len(Z0)):
                # 上下边
                if i <=1 or i >=6:#内表面边，上下边
                    number = int(p.edges[Z0[i][1]].getSize(printResults=False)/3)  # 获取边长度
                    if number == 0 :
                        number = 1
                    upANDdownEdges = partedges.findAt((Z0[i][0],))
                    p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 左右边布种子
                # 左右边
                else:
                    leftANDraghtEdges = partedges.findAt((Z0[i][0],))
                    p.seedEdgeByNumber(edges=leftANDraghtEdges, number=1)  # 左右边布种子
        elif tsX[0] > Xx0 or Xx0> tsX[1]:  # 封头
            for i in range(len(Z0)):
                partedges = p.edges
                try:  # 先尝试曲边寻找
                    p.edges[Z0[i][1]].getCurvature(0.5)
                    cellshangxiazong.append((Z0[i][0],))  # 曲边即是上下
                    number = int(p.edges[Z0[i][1]].getSize(printResults=False))  # 获取边长度
                    if number == 0 :  # 下面是对边种子点数量进行处理
                        number = 1
                    if number > 0:
                        number = int(number)
                    upANDdownEdges = partedges.findAt((Z0[i][0],))
                    p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 上下边布种子
                except:
                    cellzuoyouzong.append((Z0[i][0],))  # 直边即是左右
                    number = int(p.edges[Z0[i][1]].getSize(printResults=False))  # 获取边长度
                    leftANDrightEdges = partedges.findAt((Z0[i][0],))
                    p.seedEdgeByNumber(edges=leftANDrightEdges, number=1)  # 左右边布种子
        # 弧边
        i=0
        if i ==0:
            for i in range(len(ZHubian)):
                cellcurveszong.append((ZHubian[i][0],))
                number = int(p.edges[Z0[i][1]].getSize(printResults=False) / 5)  # 获取边长度并除以3得到单元布置数量初值
                if number == 0:
                    number = 1
                CurvesEdges = partedges.findAt((Z0[i][0],))
                p.seedEdgeByNumber(edges=CurvesEdges, number=number)  # 左右边布种子
            i += 1

    # 网格绘制
    p.generateMesh()

    # 为各边不同的网格控制方式生成相应的set集
    # 厚度方向网格数量控制
    partedges = p.edges
    leftANDrightEdges = partedges.findAt(cellzuoyouzong[0])
    for ii in range(1,len(cellzuoyouzong)):
        Edge= partedges.findAt(cellzuoyouzong[ii])
        leftANDrightEdges =leftANDrightEdges+Edge
    # 创建左右边集合
    p.Set(name='Set'+'-'+'HouduEdges',edges=leftANDrightEdges)
    # p.seedEdgeByNumber(edges=leftANDrightEdges,number=1,constraint=FINER)
    # 经线方向上网格数量控制
    upANDdownEdges = partedges.findAt(cellshangxiazong[0])
    for ii in range(1,len(cellshangxiazong)):
        Edge= partedges.findAt(cellshangxiazong[ii])
        upANDdownEdges =upANDdownEdges+Edge
    # 创建上下边集合
    p.Set(name='Set'+'-'+'JingxianEdges',edges=upANDdownEdges)
    # p.seedEdgeByNumber(edges=upANDdownEdges,number=5,constraint=FINER)#非筒身布种子
    # 旋转弧线上的网格数量控制
    CurvesEdges = partedges.findAt(cellcurveszong[0])
    for ii in range(1,len(cellcurveszong)):
        Edge= partedges.findAt(cellcurveszong[ii])
        CurvesEdges =CurvesEdges+Edge
    # 创建弧线边集合
    p.Set(name='Set'+'-'+'HuxianEdges',edges=CurvesEdges)
    # p.seedEdgeByNumber(edges=CurvesEdges,number=int(zhuanjiao/2),constraint=FINER)

    # 单元类型选择
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
        kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
        hourglassControl=DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)

    # 选择cell
    c = p.cells
    Regionszhiding=c.findAt((face_coordinate[0],))
    for k in range(1,len(face_coordinate)):#因为前面只用到了z=0的侧面，所以可以通过该面的中点得到cells
        Regionszhiding1=c.findAt((face_coordinate[k],))
        Regionszhiding=Regionszhiding+Regionszhiding1
    cell=p.Set(name='Set' + '-' + 'cell', cells=Regionszhiding)
    p.setElementType(regions=cell, elemTypes=(elemType1, elemType2,elemType3))

    # 转换叠层上下表面
    f = p.faces
    p.assignStackDirection(referenceRegion=f.findAt(coordinates=StackDirectionSurface), cells=Regionszhiding)


    ######################################################################################################################


    # 切换到JOB模块
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    session.viewports['Viewport: 1'].view.fitView()

    # 创建一个job
    mdb.Job(name='Job-1-composite-PressureVessel', model='composite',
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)

    # 提交job
    # mdb.jobs['Job-1-composite-PressureVessel'].submit(consistencyChecking=OFF)


    ######################################################################################################################


    # # # 转入后处理模块

     
 
 
 
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
     
# 下面是一体化的
######################################################################################################################



# 自定义数据处理相关函数
#声明自定义函数
def yitihuaChaiFen(data):
    a = data
    Line = zeros((1, 4))
    Aa = zeros((len(a), 4))
    for i in range(len(a)):  # 循环赋值
        Line0 = a[i].split()
        Aa[i][0] = float(Line0[0])
        Aa[i][1] = float(Line0[1])/2#将文件中的直径值转换为半径值
        Aa[i][2] = float(Line0[2])
        Aa[i][3] = float(Line0[3])
    Aa = Aa[Aa[:, 0].argsort()]  # 按照第1列对行排序
       # Rmax = max(Aa[0:, 1])  # 获取最大值
    # tmin = min(Aa[0:, 3])  # 获取最小值
    # c = np.where(Aa == max(Aa[0:, 1]))  # 获取直径最大值所在行号和列号，即是筒身段的直径值
    # return [Aa,c,Rmax,tmin]
    return [Aa]

# def yitihuaDataSolve(data,changdu,fengtoucha):
def yitihuaDataSolve(data):
    chaifen0=yitihuaChaiFen(data)#拆分出各个值
    Aa=chaifen0[0]

#预定义
    # 轮廓绘制用
    X =[]
    R =[]
    alpha =[]
    h = []
    # 角度和厚度用
    X11 = []
    R11 = []
    alpha11 = []
    h11 = []

#循环赋值
    for i1 in range(len(Aa)):
        # 画轮廓
        X.append(Aa[i1][0]);
        R.append(Aa[i1][1]);#半径值，文件里面是直径值，已在拆分函数里面转换为半径
        # alpha.append(Aa[i1][2]);
        # h.append(Aa[i1][3]);

        # 赋角度和厚度
        # if Aa[i1][2]>0:#只提取一个循环的数据
        #     X11.append(Aa[i1][0]);
        #     R11.append(Aa[i1][1]);  # 半径值，文件里面是直径值，已在拆分函数里面转换为半径
        #     alpha11.append(Aa[i1][2]);
        #     h11.append(Aa[i1][3]);

        # 输出全部厚度
        X11.append(Aa[i1][0]);
        R11.append(Aa[i1][1]);  # 半径值，文件里面是直径值，已在拆分函数里面转换为半径
        alpha11.append(Aa[i1][2]);
        h11.append(Aa[i1][3]);

# 获取重复值行号
    hanghao= []
    for ii1 in range(len(Aa) - 1):
        if Aa[ii1 + 1][0] == Aa[ii1][0] and abs(Aa[ii1 + 1][2]) == abs(Aa[ii1][2]):  # 位同角同(取绝对值)
            hanghao.append(ii1)
        else:
            pass
        if Aa[ii1 + 1][0] == Aa[ii1][0] and abs(Aa[ii1 + 1][2]) != abs(Aa[ii1][2]):  # 位同角不同(取绝对值)
            hanghao.append(ii1)
        else:
            pass

    hanghao1 = []
    for ii1 in range(len(X11) - 1):
        if X11[ii1 + 1] == X11[ii1] and abs(alpha11[ii1 + 1]) == abs(alpha11[ii1]):  # 位同角同(取绝对值)
            hanghao1.append(ii1)
        else:
            pass
        if X11[ii1 + 1] == X11[ii1] and abs(alpha11[ii1 + 1]) != abs(alpha11[ii1]):  # 位同角不同(取绝对值)
            hanghao1.append(ii1)
        else:
            pass

# 去除重复坐标值
    i=0
    for num in hanghao:
        num=num-i
        i=i+1
        del (X[num])
        del (R[num])
        # del (alpha[num])
        # del (h[num])
    i = 0
    for num in hanghao1:
        num = num - i
        i = i + 1
        del (X11[num])
        del (R11[num])
        del (alpha11[num])
        del (h11[num])

# 2020.8.23改，用总长来限定这里的尾部下降段
#     weizhi=[i for i in range(len(X)) if X[i]>(changdu-fengtoucha) and (R[i]-R[i-1])<0]
#     X=X[0:weizhi[0]]
#     R=R[0:weizhi[0]]


    weizhi=np.where(X11<=max(X))
    # weizhi=max([i for i in range(len(X11)) if X11[i]<=X[weizhi[-1]]])
    X11 = X11[0:weizhi[0][-1]+1]
    R11 = R11[0:weizhi[0][-1]+1]
    alpha11 = alpha11[0:weizhi[0][-1]+ 1]
    h11 = h11[0:weizhi[0][-1] + 1]

    return [X,R,X11,R11,alpha11,h11]###X11,alpha11原始LAM文件位置角度


 #BSpline3次
# degree. 3次样条
    # 第一步：设置曲线次数为p = 3，读取控制顶点P0, P1, ...Pn，根据控制顶点得到节点向量U = {u0, u1....um}，注意m = n + p + 1.
    # 为了保证拟合的曲线经过第一个、最后一个控制点，节点向量首末重复度要设置为为p + 1，即U = {0, 0, 0, 0, u(4)...u(n), 1, 1, 1, 1}
    # 第二步：按照本文描述的DeBoor算法完成递推函数
    # 第三步：在0 ≤ u ≤ 1取值，即可得到u对应的P(u)位置
# 规范化初始弧长
def yitihuaguifanhua(CPointX, CPointY):
    knotType = 1  # knot的计算方式 0：平均取值  1：根据长度取值
    k = 3
    n = len(CPointX)
    knot = list(nmp.arange(n + k + 1))  # 创建首末控制点重度，便于曲线过端点控制点
    for i in range(0, 4):  # 由B样条性质构建左右数据控制点的k+1重节点，k为基函数阶数
        knot[i] = 0.0
        knot[n + i] = 1.0
    if knotType:
        L = list(nmp.arange(n - 1))
        S = 0
        for i in range(n - 1):  # 计算矢量长
            L[i] = nmp.sqrt(pow(CPointX[i + 1] - CPointX[i], 2) + pow(CPointY[i + 1] - CPointY[i], 2))
            S = S + L[i]
        tmp = L[0]
        for i in range(4, n):  # 按矢量比例均分，曲线定义域内
            tmp = tmp + L[i - 3]
            knot[i] = tmp / S
    else:
        tmp = 1 / (float(n) - 3)  # 因为三次样条是4个节点一个曲线段，现在在将曲线定义域用n-k来均分定义区间
        tmpS = tmp
        for i in range(4, n):
            knot[i] = tmpS
            tmpS = tmpS + tmp
    return knot
# 基函数值递归计算
def yitihuacoxDeBoor(u, knots, i, k):
    # Test for end conditions
    if (k == 0):
        if (knots[i] <= u and u < knots[i + 1]):
            return 1
        return 0
    Den1 = knots[i + k] - knots[i]
    Den2 = knots[i + k + 1] - knots[i + 1]
    Eq1 = 0;
    Eq2 = 0;
    if Den1 > 0:
        Eq1 = ((u - knots[i]) / Den1) * yitihuacoxDeBoor(u, knots, i, (k - 1))
    if Den2 > 0:
        Eq2 = ((knots[i + k + 1] - u) / Den2) * yitihuacoxDeBoor(u, knots, (i + 1), (k - 1))

    return Eq1 + Eq2
# 生成Bspline点坐标，并修正R值
def yitihuaB3Spline(CPointX,CPointY,knot,Bsjindu):
    BsplineX = []
    BsplineY = []
    for ip in range(len(CPointX) - 3):  # 创建B样条
        u = linspace(knot[ip + 3], knot[ip + 4], Bsjindu)  # 节点区间值
        for j in arange(len(u)):  # 创建均B样条
            BsplineX.append(
                CPointX[ip] * yitihuacoxDeBoor(u[j], knot, ip, 3) + CPointX[ip + 1] * yitihuacoxDeBoor(u[j], knot, ip + 1, 3) \
                + CPointX[ip + 2] * yitihuacoxDeBoor(u[j], knot, ip + 2, 3) + CPointX[ip + 3] * yitihuacoxDeBoor(u[j], knot, ip + 3,
                                                                                                   3))
            BsplineY.append(
                CPointY[ip] * yitihuacoxDeBoor(u[j], knot, ip, 3) + CPointY[ip + 1] * yitihuacoxDeBoor(u[j], knot, ip + 1, 3) \
                + CPointY[ip + 2] * yitihuacoxDeBoor(u[j], knot, ip + 2, 3) + CPointY[ip + 3] * yitihuacoxDeBoor(u[j], knot, ip + 3,
                                                                                                   3))
    # 2020.8.30改
    # 会导致出现点过原点的错误
    if BsplineX[-1] == 0:
        del BsplineX[-1]
        del BsplineY[-1]
    # 到尾喷管端点
    BsplineX.append(CPointX[-1])
    BsplineY.append(CPointY[-1])
    # plt.plot(BsplineX, BsplineY, '-r')
    # plt.plot([BsplineX[0], BsplineX[-1]], [0, 0], 'g')
    # plt.plot([BsplineX[0], BsplineX[-1]], [R0_jikong, R0_jikong], 'b')
    # plt.plot(CPointX, CPointY, '-ob', label="Waypoints")
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()
    return [BsplineX,BsplineY]

#通过厚度求取外轮廓，即是通过两个数据点的向量值夹角和厚度值反求对应的外轮廓点
def yitihuawailunkuo1(X,R,H):
    X2wai = []
    R2wai = []

# 第一个数据点位的厚度反算轮廓
    theta = math.atan(-(X[1] - X[0]) / (R[1] - R[0]))
    Delta_X = (2 * H[0]) * cos(theta)
    Delta_R = (2 * H[0]) * sin(theta)
    X2wai.append(X[0] - Delta_X)
    R2wai.append(R[0] - Delta_R)

# 除第一个数据点外的数据点位厚度反算轮廓
#     2020.10.1改
    for i111 in range(len(X) - 1):
        i111 = i111 + 1
        theta=math.atan(-(R[i111]-R[i111-1])/(X[i111]-X[i111-1]))
        Delta_X=(2*H[i111])*sin(theta)
        Delta_R=(2*H[i111])*cos(theta)
        if R[i111]>R[i111-1]:
            X2wai.append(X[i111]+Delta_X)
            R2wai.append(R[i111]+Delta_R)
        if R[i111]<R[i111-1]:
            X2wai.append(X[i111]+Delta_X)
            R2wai.append(R[i111]+Delta_R)
        if R[i111]==R[i111-1]:
            X2wai.append(X[i111])
            R2wai.append(R[i111]+ 2*H[i111])
    # plt.plot(X2wai,R2wai,'g')
    # plt.show()
    return [X2wai,R2wai]

#通过厚度求取外轮廓，即是通过两个数据点的向量值夹角和厚度值反求对应的外轮廓点
#左封头
def yitihuawailunkuozuo(X,R,H,tsX):
    X2_zuo_wai = list()
    R2_zuo_wai = list()
    # if X[0]!=max(X):
    #     X=X[::-1]#因为以赤道圆作为初始位置，需翻转左封头的数据
    #     R=R[::-1]
    #     H=H[::-1]
    X2_zuo_wai.append(tsX[0])
    R2_zuo_wai.append(max(R)+2*H[0])
    for i111 in range(len(X) - 1):
        i111 = i111 + 1
        theta=math.atan(-(X[i111]-X[i111-1])/(R[i111]-R[i111-1]))
        Delta_X=(2*H[i111])*cos(theta)
        Delta_R=(2*H[i111])*sin(theta)
        X2_zuo_wai.append(X[i111]-Delta_X)
        R2_zuo_wai.append(R[i111]-Delta_R)
    # plt.plot(X2_zuo_wai,R2_zuo_wai)
    return [X2_zuo_wai,R2_zuo_wai]

#右封头
def yitihuawailunkuoyou(X,R,H,tsX):
    X2_you_wai = list()
    R2_you_wai = list()

    X2_you_wai.append(tsX[1])
    R2_you_wai.append(max(R)+ 2*H[0])
    for i222 in range(len(X) - 1):
        i222 = i222 + 1

        theta = math.atan(-(R[i222] - R[i222 - 1]) / (X[i222] - X[i222 - 1]))
        Delta_X = (2 * H[i222]) * sin(theta)
        Delta_R = (2 * H[i222]) * cos(theta)
        if R[i222] > R[i222 - 1]:
            X2_you_wai.append(X[i222] + Delta_X)
            R2_you_wai.append(R[i222] + Delta_R)
        if R[i222] < R[i222 - 1]:
            X2_you_wai.append(X[i222] + Delta_X)
            R2_you_wai.append(R[i222] + Delta_R)
        if R[i222] == R[i222 - 1]:
            X2_you_wai.append(X[i222])
            R2_you_wai.append(R[i222] + 2 * H[i222])
            R2_you_wai.append(R[i222] + 2 * H[i222])

#         k1=-(X[i222] - X[i222 - 1]) / (R[i222] - R[i222 - 1])
#         theta=math.atan(k1)
#         Delta_X=(2*H[i222])*cos(theta)
#         if k1>=0:
#             Delta_R=(2*H[i222])*sin(theta)
#         if k1<0:#喉颈部位以后的厚度反推计算
#             Delta_R = -(2 * H[i222]) * sin(theta)
# # 一体化独有
#         if R[i222]>R[i222-1]:
#             X2_you_wai.append(X[i222]-Delta_X)
#             R2_you_wai.append(R[i222]+Delta_R)
#         if R[i222]<R[i222-1]:
#             X2_you_wai.append(X[i222]+Delta_X)
#             R2_you_wai.append(R[i222]+Delta_R)
#         if R[i222]==R[i222-1]:
#             X2_you_wai.append(X[i222])
#             R2_you_wai.append(R[i222]+ 2*H[i222])
    # plt.plot(X2_you_wai,R2_you_wai,'r')
    # plt.show()
    return [X2_you_wai,R2_you_wai]

def yitihuaXiuZheng(x,X,R):#R坐标修正
    # 修正半径值,因为X/Y坐标是一一对应的，所以可以找到离x坐标最近的那个R作为我的现在的修正的R,为后面厚度反推外轮廓做准备
    Rxiuzh = []
    Xxiuzh = []
    Hxiuzh = []
    for i in range(len(x)):
        jdz = [abs(X[j] - x[i]) for j in range(len(X))]  # 绝对值
        c1 = jdz.index(min(jdz))  # 获取各个x对应为修正位置
        Xxiuzh.append(X[c1])
        Rxiuzh.append(R[c1])
    return [Xxiuzh,Rxiuzh]

def yitihuahoudupingwen(X,R,pingwenxishu,tsX):#前后前跳动值平稳处理，用来为轮廓生成做准备
    # 两个向量夹角判定
    j=0
    weizhi=[]
    for i in range(len(X) - 2):
        axaingliang=[X[i+1-j]-X[i-j],R[i+1-j]-R[i-j]]#a向量
        bxaingliang=[[X[i+2-j]-X[i+1-j]],[R[i+2-j]-R[i+1-j]]]#b向量
        amochang=((X[i+1-j]-X[i-j])**2+(R[i+1-j]-R[i-j])**2)**0.5#模长
        bmochang=((X[i+2-j]-X[i+1-j])**2+(R[i+2-j]-R[i+1-j])**2)**0.5#模长
        alp=(math.acos((np.dot(axaingliang,bxaingliang)) / (amochang* bmochang)))*180/math.pi
        if abs(alp)>pingwenxishu*180/math.pi and X[i-j]>tsX[1]:#去除急转的角度相邻位置值
            weizhi.append(i+2)#将波动位置存储起来后面统一删除
        if abs(alp)>pingwenxishu*180/math.pi and X[i-j]<tsX[0]:#去除急转的角度相邻位置值
            weizhi.append(i+2)#将波动位置存储起来后面统一删除
    for ii in range(len(weizhi)):
        del (X[weizhi[ii]-j])
        del (R[weizhi[ii]-j])
        j+=1

    return [X,R]

def yitihuaQujiange(old_listX,old_listR,shuliang):
    new_listX = []
    new_listR = []
    new_listX.append(old_listX[0])
    new_listR.append(old_listR[0])
    # 获取平均距离
    distance = [((old_listX[i] - old_listX[i + 1]) ** 2 + (old_listR[i] - old_listR[i + 1]) ** 2) ** 0.5 for i in
               range(len(old_listX) - 1)]
    pingjudistance = sum(distance) / shuliang
    i = 0
    for j in range(1, len(old_listX)):
        if ((old_listX[j] - old_listX[i]) ** 2 + (old_listR[j] - old_listR[i]) ** 2) ** 0.5 > pingjudistance:
            new_listX.append(old_listX[j - 1])
            new_listR.append(old_listR[j - 1])
            i = j
    # new_list=old_list[0:len(old_list):inter]
    if new_listX[-1] != old_listX[-1] or new_listR[-1] != old_listR[-1]:
        new_listX.append(old_listX[-1])
        new_listR.append(old_listR[-1])
    else:
        pass
    return [new_listX, new_listR]

# 2020.8.23修正增加，这里的x会出现相等的情况必须再次去重处理
def yitihuaquchong(X,R,H,alpha1):
    weizhi=[i for i in range(len(X)) if X[0:i+1].count(X[i])>1]# 找到去除X重复项的位置
    for i in range(len(weizhi)):
        del (X[weizhi[i]-i])
        del (R[weizhi[i]-i])
        del (H[weizhi[i]-i])
        del (alpha1[weizhi[i]-i])
    return [X,R,H,alpha1]

# 厚度光滑性处理
def yitihuasmooth(H):
    x=arange(0,len(H),1)
    # plt.plot(x,H,'-r')
    H1=H[:]
    # 2020.11.6先来一个简单粗暴的处理，取两点平均
    for i in arange(1,len(H1)-1,1):
        H1[i]=(H1[i-1]+H1[i+1])/2
    H1[i+1]=(H1[i-1]+H1[i])/2
    H1[0]=(H1[0]+H1[1])/2
    # 得到内部跳跃的最大差值
    maxtiaoyuechazhi0=max([H1[i+1]-H1[i] for i in range(len(H1)-1)])
    a=20
    tiaoyue=[]
    tiaoyue.append(maxtiaoyuechazhi0)
    houducha=[]
    xielv0 = []
    for i in arange(0, a, 1):
        for i in arange(1, len(H1) - 1, 1):
            H1[i] = (H1[i - 1] + H1[i + 1]) / 2
        H1[i + 1] = (H1[i - 1] + H1[i]) / 2
        H1[0] = (H1[0] + H1[1]) / 2
        maxtiaoyuechazhi1 = max([H1[i + 1] - H1[i] for i in range(len(H1) - 1)])
        tiaoyue.append(maxtiaoyuechazhi1)
        chazhi = [abs(H1[i] - H[i]) for i in range(len(H1))]
        zuidazhi = max(chazhi)
        houducha.append(zuidazhi)
        xyz = [tiaoyue[i + 1] / houducha[i] for i in range(len(tiaoyue) - 1)]
        if i > 0:
            xielv = (xyz[i] - xyz[i - 1]) / 1#通过斜率控制，刚开始下降比较快，后面就没必要了，误差还大
            xielv0.append(xielv)
            if abs(xielv) < 2e-3:
                break
    # plt.plot(arange(0, len(H1), 1),H1, '-b')
    # # plt.plot(arange(0, len(tiaoyue), 1), tiaoyue, '-b')
    # # plt.plot(arange(1, len(houducha) + 1, 1), houducha, '-r')
    # # plt.plot(arange(1, len(xyz) + 1, 1), xyz, '-g')
    # plt.grid(True)
    # plt.legend()
    # # plt.axis("equal")
    # plt.show()

    return H1






######################################################################################################################

# def yitihua(fileName,layup,E1,E2,E3,G12,G13,G23,Nu12,Nu13,Nu23,pressure,hxdt,
        # R0_jikong,R0_weipenguan,changdu,fengtoucha,youfengtou,tschangdu,zhuanjiao,jinshunei):
def yitihua(fileName,layup,E1,E2,E3,G12,G13,G23,Nu12,Nu13,Nu23,pressure,hxdt,changdu,tschangdu,zhuanjiao,jinshunei):
    print(layup)
    # print(E1,E2,E3,G12,G13,G23,Nu12,Nu13,Nu23,pressure,hxdt,\
        # R0_jikong,R0_weipenguan,changdu,fengtoucha,youfengtou,tschangdu,zhuanjiao,jinshunei)
    print(E1,E2,E3,G12,G13,G23,Nu12,Nu13,Nu23,pressure,hxdt,changdu,tschangdu,zhuanjiao,jinshunei)
    
    # 导入相关的cadwind数据
    # 读取数据文件
    Lam1 = open(fileName, 'r')
    data1 = Lam1.readlines()
    Lam1.close()
    del (data1[0:2])  # 删除前面的两行非数值项,Cadwind默认生成的非数据行


    ######################################################################################################################


    puceng=[]

    for x in layup:
        puceng.append(x[0])

        

    engineeringdata=(E1,E2,E3,Nu12,Nu13,Nu23,G12,G13,G23)  

    #默认参数不需要插件输入
    #Bsjindu=BsplinePrecision#B样条单个区间拟合精度   
    # todegree=pi/180.0
    # jiangeshu=internalNum#绘制三维模型间隔取点数量
    # pingwenxishu=smoothCoefficient#外轮廓相邻点平稳系数，越小越平稳，太小会出问题，需适当调节

    Bsjindu=20#B样条单个区间拟合精度
    zhuanjiao=30.0#三维实体模型的转角
    todegree=pi/180.0
    pingwenxishu=0.15#外轮廓相邻点平稳系数，越小越平稳，太小会出问题，需适当调节


    ######################################################################################################################


    # 初始化数据存储参数
    # 端点值
    # 左
    XddianZ=[]
    RddianZ=[]
    # 右
    XddianY=[]
    RddianY=[]

    # 轮廓
    # 内表面
    neisX=[]
    neisR=[]
    # 外表面
    waisX=[]
    waisR=[]

    ######################################################################################################################


     # 得到传入处理后的初步数据
    #内表面
    # 得到最内层的轮廓和螺旋缠绕的角度数据
    # 得到各个数据值
    # data=yitihuaDataSolve(data1,changdu,fengtouchazuo)
    data=yitihuaDataSolve(data1)
    #绘制轮廓用
    X1=data[0]
    R1=data[1]
    # 角度赋予用
    X12=data[2]
    R12=data[3]
    alpha1=data[4]
    h1=data[5]


    ######################################################################################################################


    # 2020.11.7
    # 极孔线性外插
    # 左
    x1=X12[0]
    x2=X12[1]

    y1=R12[0]
    y2=R12[1]

    Y=np.array([[y1],[y2]])
    mat2=np.array([[x1,1],[x2,1]])
    mat2 = np.mat(mat2)
    mat3 = mat2.I#求逆
    xishu=mat3*Y
    zuojikquzhengX=-ceil(abs(X12[0]))
    zuojikquzhengR=round(xishu[0,0]*zuojikquzhengX+xishu[1,0])

    R0_jikong=zuojikquzhengR
    fengtouchazuo=abs(zuojikquzhengX)#左封头长度

    # 去除尾部下降段
    j=0
    for i in range(len(X1)):
        if X1[i-j]>changdu-fengtouchazuo:
            del (X1[i-j])
            del (R1[i-j])
            del (X12[i-j])
            del (R12[i-j])
            del (alpha1[i-j])
            del (h1[i-j])
            j+=1

    ######################################################################################################################


    # 根据两边极孔条件删除不在范围内的点
    weizhi00=[i for i in range(len(R12)) if R12[i]<=R0_jikong and X12[i]<0]
    j=0
    for i in range(len(weizhi00)):
        del (X12[weizhi00[i-j]])
        del (R12[weizhi00[i-j]])
        del (alpha1[weizhi00[i-j]])
        del (h1[weizhi00[i-j]])
        j+=1
    aa=[i for i in range(len(X12)) if X12[i]<=0]
    if aa[-1]!=0:
        X12[aa[-1]]=0
        R12[aa[-1]]=R12[aa[-1]+1]

    ######################################################################################################################


    # # 对厚度进行平稳处理，便于后面轮廓的光滑
    h1=yitihuasmooth(h1)

    CPointX=[]
    CPointY=[]
    # 用B样条拟合外轮廓
    CPointX = X1[:]  # 整理B样条需要数据
    CPointY = R1[:]
    # 根据两边极孔条件删除不在范围内的点
    weizhi00=[i for i in range(len(CPointY)) if CPointY[i]<=R0_jikong and CPointX[i]<0]
    j=0
    for i in range(len(weizhi00)):
        del (CPointX[weizhi00[i-j]])
        del (CPointY[weizhi00[i-j]])
        j+=1

    ######################################################################################################################


    # 右侧插入控制点减缓曲线突变程度
    p=20#插入控制点数量
    qujian=20
    for i in range(qujian):
        for j in arange(1,p,1):
            [CPointX.insert(-(p*i+j),X1[-(i+1)]-j*(X1[-(i+1)]-X1[-(i+2)])/p)]
            [CPointY.insert(-(p*i+j),R1[-(i+1)]-j*(R1[-(i+1)]-R1[-(i+2)])/p)]


    ######################################################################################################################


    # 得到B样条的数据点
    knot = yitihuaguifanhua(CPointX[:], CPointY[:])  # 生成B样条节点
    B3 = yitihuaB3Spline(CPointX[:], CPointY[:], knot[:],Bsjindu)  # 得到B样条数据点

    # 2020.9.19改
    BsplineX = B3[0][:]
    BsplineY = B3[1][:]
    x1 = BsplineX[0] + (R0_jikong - BsplineY[1]) * ((BsplineX[0] - BsplineX[1]) / (BsplineY[0] - BsplineY[1]))
    BsplineX.insert(0, x1)
    BsplineY.insert(0, R0_jikong)
    B3[0][:] = BsplineX
    B3[1][:] = BsplineY


    ######################################################################################################################

    # 内轮廓
    neisX=B3[0][:]
    neisR=B3[1][:]
    # plt.plot(neisX,neisR)
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()
    # 端点
    # 左
    XddianZ.append(B3[0][0])
    RddianZ.append(B3[1][0])
    # 右
    XddianY.append(B3[0][-1])
    RddianY.append(B3[1][-1])


    ######################################################################################################################


    #外表面
    # 得到螺旋缠绕的外表面轮廓数据
    # 修正现在的X,R值，用B样条中的离该点X坐标最近的那个数据点代替原始X,R值
    cc = yitihuaXiuZheng(X12, B3[0], B3[1])
    xiuzX = cc[0][:]
    xiuzR = cc[1][:]


    bb=yitihuaquchong(xiuzX[:], xiuzR[:], h1[:],alpha1[:])
    xiuzX =bb[0][:]
    xiuzR = bb[1][:]
    h1= bb[2][:]
    alpha1= bb[3][:]


    ######################################################################################################################


    # 得到环向层的厚度之后加入到螺旋层的筒身厚度计算中去
    huanxsl=puceng.count('Hoop')
    tshouduhx=huanxsl*hxdt
    # 螺旋总厚度
    luoxuansl=puceng.count('HelixLoop')
    h2=[i*luoxuansl for i in h1]


    ######################################################################################################################


    # 后面需要处理的厚度
    XX11=[]
    RR11=[]
    # 找出筒身的位置
    tsX=[]
    tsR=[]
    XX =fengtouchazuo  ###左边封头位置
    zuotongshen=0
    jdz = [abs(B3[0][j] - XX) for j in range(len(B3[0]))]  # 绝对值
    weizhi1=jdz.index(min(jdz))  # 获取各个x对应为修正位置
    #得到左部赤道圆位置坐标
    # tsX.append(B3[0][weizhi1])
    # tsR.append(B3[1][weizhi1])
    tsX.append(0)
    tsR.append(B3[1][weizhi1])

    # 提取筒身螺旋单层厚度/角度
    weizhi=[i for i in range(len(xiuzX)) if 0<xiuzX[i]<tschangdu]
    lxdt=h1[(weizhi[0]+weizhi[-1])/2]
    luoxuanjiao=abs(alpha1[(weizhi[0]+weizhi[-1])/2])

    # 得到右部赤道圆位置坐标
    XX =tsX[0] + tschangdu     #####右边赤道圆位置
    youtongshen=tsX[0]
    jdz = [abs(B3[0][j] - XX) for j in range(len(B3[0]))]  # 绝对值
    weizhi2=jdz.index(min(jdz))  # 获取各个x对应为修正位置
    # tsX.append(B3[0][weizhi2])
    tsR.append(B3[1][weizhi2])
    tsX.append(tsX[0] + tschangdu)
    # tsR.append(tsR[0])#保证两侧赤道圆直径相等

    # 左
    a=[xiuzX[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0] and xiuzX[i]<0]
    b=[xiuzR[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0] and xiuzX[i]<0]
    c=[h2[i] for i in range(len(xiuzX)) if xiuzX[i]<tsX[0] and xiuzX[i]<0]
    X=list(reversed(a))
    R=list(reversed(b))
    H1=list(reversed(c))
    waizuo=yitihuawailunkuozuo(X,R,H1,tsX)#得到极孔处修正的外轮廓，删掉了一部分厚度
    Xzuo=list(reversed(waizuo[0]))
    Rzuo=list(reversed(waizuo[1]))
    Hzuo=list(reversed(H1))
    j=0
    for i in range(len(Xzuo)):
        if Xzuo[i-j] <X[0]-fengtouchazuo-lxdt*luoxuansl*max(Rzuo)/R0_jikong:#长度限定,对内外轮廓多做截断处理
            #外部
            del(Xzuo[i-j])
            del(Rzuo[i-j])
            del(Hzuo[i-j])
            del(a[i-j])
            del(b[i-j])
            j+=1
        if Xzuo[i-j]>tsX[0]:
            del (Xzuo[i - j])
            del (Rzuo[i - j])
            del (Hzuo[i - j])
            del (a[i - j])
            del (b[i - j])
            j += 1
    [XX11.append(a[i]) for i in range(len(a)) if (a[i]<=max(Xzuo)) and (a[i]>=min(Xzuo))]
    [RR11.append(b[i]) for i in range(len(a)) if (a[i]<=max(Xzuo)) and (a[i]>=min(Xzuo))]
    ############################################################
    # 得到右部赤道圆位置坐标
    XX =tsX[0] + tschangdu     #####右边赤道圆位置
    youtongshen=XX
    jdz = [abs(B3[0][j] - XX) for j in range(len(B3[0]))]  # 绝对值
    weizhi2=jdz.index(min(jdz))  # 获取各个x对应为修正位置
    # tsX.append(B3[0][weizhi2])
    # tsR.append(B3[1][weizhi2])

    # 右
    X=[xiuzX[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    R=[xiuzR[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    Hyou=[h2[i] for i in range(len(xiuzX)) if xiuzX[i]>tsX[1]]
    waiyou=yitihuawailunkuoyou(X,R,Hyou,tsX)#得到极孔处修正的外轮廓，删掉了一部分厚度
    Xyou=waiyou[0]
    Ryou=waiyou[1]
    j=0
    for i in range(len(Xyou)):
        if Xyou[i-j] >changdu-abs(min(Xzuo)):#长度限定
           # 外部
            del (Xyou[i - j])
            del (Ryou[i - j])
            del (Hyou[i - j])
            del (X[i - j])
            del (R[i - j])
            j += 1
        if Xyou[i-j]<tsX[1]:
            del (Xyou[i - j])
            del (Ryou[i - j])
            del (Hyou[i - j])
            del (X[i - j])
            del (R[i - j])
            j += 1
    # 筒身数据
    [XX11.append(xiuzX[i]) for i in range(len(xiuzX)) if (xiuzX[i]>=tsX[0]) and (xiuzX[i]<=tsX[1])]
    [RR11.append(xiuzR[i]) for i in range(len(xiuzX)) if (xiuzX[i]>=tsX[0]) and (xiuzX[i]<=tsX[1])]
    # 右封头
    [XX11.append(X[i]) for i in range(len(X)) if (X[i]>=min(Xyou)) and (X[i]<=max(Xyou))]
    [RR11.append(R[i]) for i in range(len(X)) if (X[i]>=min(Xyou)) and (X[i]<=max(Xyou))]


    ######################################################################################################################


    # 删除掉极孔处有问题的数据点之后，重新组装新的外部轮廓数据
    #重新定义新的厚度数据，保护原始数据
    X11=[]
    R11=[]
    H11=[]

    #左封头数据
    [X11.append(i) for i in Xzuo]
    [R11.append(i) for i in Rzuo]
    [H11.append(i) for i in Hzuo]

    #筒身数据
    for i in range(len(xiuzX)):
        if xiuzX[i]<=tsX[1] and xiuzX[i]>=tsX[0]:
            X11.append(xiuzX[i])
            R11.append(xiuzR[i])#一体化的半径有的时候两侧半径大小不一问题，后面可能会出错
            # R11.append(tsR[0]+tshouduhx)
            # H11.append(h2[i])
            H11.append(lxdt*luoxuansl)

    #右封头数据
    [X11.append(i) for i in Xyou]
    [R11.append(i) for i in Ryou]
    [H11.append(i) for i in Hyou]


    ######################################################################################################################


    # 做赤道圆的过渡处理
    # 2020.9.22
    # 在赤道圆过渡位置寻找使其与封头相切的位置近似点
    leftTransitionX=tsX[0]
    leftTransitionR=tsR[0]+lxdt*luoxuansl*2+tshouduhx
    rightTransitionX=tsX[1]
    rightTransitionR=tsR[0]+lxdt*luoxuansl*2+tshouduhx
    tangentlineL=[]
    tangentlineR=[]
    wz=[]
    for i in range(len(X11)):
        if X11[i]<leftTransitionX:
            tangentlineL.append(math.atan(((X11[i] -leftTransitionX) / (R11[i] -leftTransitionR)))*180/math.pi)
        if X11[i] > rightTransitionX:
            tangentlineR.append(math.atan(((X11[i] -rightTransitionX)/(R11[i] - rightTransitionR)))*180/math.pi)
        if (X11[i] <= rightTransitionX) and (X11[i]>=leftTransitionX):
            wz.append(i)
    wz1=tangentlineL.index(max(tangentlineL))
    # wz2=tangentlineR.index(min(tangentlineR))
    # 2020.11.7
    pp=[]
    for i in range(len(tangentlineR)-1):
        xielv=tangentlineR[i+1]-tangentlineR[i]
        pp.append(xielv)
        if xielv>0:
            break
    wz2=i
    # 左
    for num in arange(wz1,min(wz)-1,1):
        t1=(tshouduhx)/(leftTransitionX-X11[wz1])*(X11[num]-X11[wz1])
        H11[num]=H11[num]+t1*0.5

    # 右
    for num in arange(max(wz)+1,max(wz)+1+wz2,1):
        t1=(tshouduhx)/(rightTransitionX-X11[max(wz)+1+wz2])*(X11[num]-X11[max(wz)+1+wz2])
        H11[num]=H11[num]+t1*0.5

    # 筒身部分环向层做加上环向层厚度的处理
    for num in arange(min(wz)-1,max(wz)+1,1):
        H11[num]=H11[num]+tshouduhx*0.5


    ######################################################################################################################


    # 对厚度光滑处理减缓厚度突变
    Hh=H11
    # Hh=yitihuasmooth(H11)#光滑处理一下
    Xnei=[]
    Rnei=[]
    for ii in range(len(XX11)):#通过内部轮廓数据截断外部X两侧，后面做修正处理
        if XX11[ii]<=max(X11) and XX11[ii]>=min(X11):
            Xnei.append(XX11[ii])
            Rnei.append(RR11[ii])
    try:
        wailunkuo = yitihuawailunkuo1(Xnei[:], Rnei[:], Hh[:])
    except:
        print("the wailunkuo command have problem")
        # 去除掉外轮廓x坐标在内轮廓横坐标控制范围外的点
    # for i in  range(len(wailunkuo[0])):
    #     plt.plot((Xnei[i],wailunkuo[0][i]),(Rnei[i],wailunkuo[1][i]),'-b')
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()

    j=0
    for ii in range(len(wailunkuo[0][:])):
        if wailunkuo[0][ii-j]<min(Xnei[:]) or  wailunkuo[0][ii-j]>max(Xnei[:]):
            del(wailunkuo[0][ii-j])
            del(wailunkuo[1][ii-j])
            j += 1
    WX = wailunkuo[0][:]
    WR = wailunkuo[1][:]

    # plt.plot(WX,WR)
    # plt.plot(Xnei,Rnei)
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()


    ######################################################################################################################


    # 厚度平稳处理
    # 左封头
    xzuo=[WX[i] for i in range(len(WX[:])) if WX[i]<=tsX[0]]
    rzuo=[WR[i] for i in range(len(WX[:])) if WX[i]<=tsX[0]]
    chuli=yitihuahoudupingwen(xzuo[::-1],rzuo[::-1],pingwenxishu,tsX)
    WXzuo = chuli[0][::-1]
    WRzuo = chuli[1][::-1]
    # 筒身
    WXzhong=[WX[i] for i in range(len(WX[:])) if tsX[0]<WX[i]<tsX[1]]
    WRzhong=[WR[i] for i in range(len(WX[:])) if tsX[0]<WX[i]<tsX[1]]
    # 右封头
    xyou=[WX[i] for i in range(len(WX[:])) if WX[i]>=tsX[1]]
    ryou=[WR[i] for i in range(len(WX[:])) if WX[i]>=tsX[1]]
    chuli=yitihuahoudupingwen(xyou[:],ryou[:],pingwenxishu,tsX)
    WX = WXzuo+WXzhong+chuli[0][:]
    WR = WRzuo+WRzhong+chuli[1][:]

    # plt.plot(WX,WR,'*r')
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()


    ######################################################################################################################


    # 这里的抛物线修正
    # 左侧
    zuoX11=min(neisX)-Hzuo[0]#左边极孔纤维堆积后的坐标
    # 最后两点的切线值
    Kk=(WX[0]-WX[1])/(WR[0]-WR[1])
    aa=(WX[0]-zuoX11-Kk*(WR[0]-R0_jikong))/(((WR[0]**2-R0_jikong**2))-2*WR[0]*(WR[0]-R0_jikong))
    bb=Kk-2*aa*WR[0]
    cc=zuoX11-aa*R0_jikong**2-bb*R0_jikong
    XXjikongzuo=[]
    RRjikongzuo=[]
    for i in arange(R0_jikong,WR[0],0.05):#靠近极孔这边的抛物线修正
        XXjikongzuo.append(aa*i**2+bb*i+cc)
        RRjikongzuo.append(i)
    WX=XXjikongzuo+WX
    WR=RRjikongzuo+WR

    # plt.plot(WX,WR)
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()


    ######################################################################################################################


    # 进行平稳处理
    pingwen=yitihuahoudupingwen(WX,WR,pingwenxishu,tsX)
    # 用B样条拟合外轮廓
    CPointX = pingwen[0][:]  # 整理B样条需要数据
    CPointY = pingwen[1][:]

    # 2020.10.1改，右侧没到最右端处理
    if CPointX[-1]<neisX[-1]:
        CPointX.append(neisX[-1])
        # CPointY.append(neisR[-1]+Hh[-1])
        # 通过斜率插入
        xielv=(CPointY[-2]-CPointY[-1])/(CPointX[-3]-CPointX[-2])
        y0=-xielv*(CPointX[-2]-CPointX[-1])+CPointY[-1]
        CPointY.append(y0)


    # 一体化独有
    # 厚度反推后存在x交错情况，可以通过判定x值调整
    j=0
    for i in range(len(CPointX)-1):
        if CPointX[i+1-j]<CPointX[i-j]:
            del(CPointX[i+1-j])
            del(CPointY[i+1-j])
            j+=1


    # # 2020.10.1改,加密控制顶点
    # 右侧插入控制点减缓曲线突变程度
    p=20#插入控制点数量
    qujian=20
    xzuo = CPointX[:]
    xyou = CPointX[:]
    yzuo = CPointY[:]
    yyou = CPointY[:]
    for i in range(qujian):

        for j in arange(1,p,1):

            [CPointX.insert(-(p*i+j),xyou[-(i+1)]-j*(xyou[-(i+1)]-xzuo[-(i+2)])/p)]
            [CPointY.insert(-(p*i+j),yyou[-(i+1)]-j*(yyou[-(i+1)]-yzuo[-(i+2)])/p)]


    ######################################################################################################################


    # 对赤道圆处的赤道数据轮廓不明显采取的赤道圆加密数据处理，得到合适赤道圆数据
    # 筒身赤道处加密数据控制点，减少这里的变化突变
    # 找到赤道圆在控制点中的位置
    where1=[i for i in range(len(CPointX)-1) if CPointX[i]<=tsX[0] and tsX[0]<=CPointX[i+1]]#左赤道圆
    where2=[i for i in range(len(CPointX)-1) if CPointX[i]<=tsX[1] and tsX[1]<=CPointX[i+1]]#右赤道圆
    CPointX.insert(where1[0]+1,tsX[0])#扩充赤道圆值
    CPointY.insert(where1[0]+1,CPointY[where1[0]+1])#扩充赤道圆值
    CPointX.insert(where2[0]+2,tsX[1])
    CPointY.insert(where2[0]+2,CPointY[where2[0]+1])#扩充赤道圆值
    where1=[i for i in range(len(CPointX)) if CPointX[i]==tsX[0]]#左赤道圆
    where2=[i for i in range(len(CPointX)) if CPointX[i]==tsX[1]]#右赤道圆


    ######################################################################################################################


    # 赤道圆过度处加密
    p=5
    qujian=[where1[0],where2[0]]
    for i in qujian:
        if i==qujian[0]:
            for j in arange(1,p,1):
                [CPointX.insert((i+1),CPointX[(i)]+j*(CPointX[(i+1)]-CPointX[(i)])/p)]
                [CPointY.insert((i+1),CPointY[(i+1)])]
        if i == qujian[1]:
            for j in arange(1, p, 1):
                [CPointX.insert((i +p-1), CPointX[(i +p-1)]-j * (CPointX[(i + p-1)] - CPointX[(i+p-2)]) / p)]
                [CPointY.insert((i +p-1), CPointY[(i - 2+p)])]

    # plt.plot(CPointX,CPointY)
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()


    ######################################################################################################################


    # 生成外轮廓B样条节点
    knot = yitihuaguifanhua(CPointX[:], CPointY[:])  # 生成B样条节点
    B3 = yitihuaB3Spline(CPointX[:], CPointY[:], knot[:],Bsjindu)  # 得到B样条数据点
    # 2020.9.19改
    BsplineX=B3[0][:]
    BsplineY=B3[1][:]
    x1 = BsplineX[0] + (R0_jikong - BsplineY[1]) * ((BsplineX[0] - BsplineX[1]) / (BsplineY[0] - BsplineY[1]))
    BsplineX.insert(0,x1)
    BsplineY.insert(0,R0_jikong)
    # 2020.10.1改，右侧没到最右端处理
    if BsplineX[-1]<neisX[-1]:
        BsplineX.append(neisX[-1])
        BsplineY.append(BsplineY[-1])


    ######################################################################################################################


    # 得到外轮廓绘制图形数据
    waisX=BsplineX
    waisR=BsplineY
    # plt.plot(waisX,waisR)
    # plt.plot(neisX,neisR)
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()
    # 端点
    # 左
    XddianZ.append(waisX[0])
    RddianZ.append(waisR[0])
    # 右
    XddianY.append(B3[0][-1])
    RddianY.append(B3[1][-1])


    # 从B样条数据里面取一部分，防止导入abaqus卡死
    neisXR=yitihuaQujiange(neisX,neisR,len(X12))
    neisX=neisXR[0][:]
    neisR=neisXR[1][:]

    waisXR=yitihuaQujiange(waisX,waisR,len(X12))
    waisX=waisXR[0][:]
    waisR=waisXR[1][:]


    # plt.plot(neisX,neisR,'*y')
    # plt.plot(waisX,waisR,'*b')
    # plt.grid(True)
    # plt.legend()
    # plt.axis("equal")
    # plt.show()
    #
    # plt.plot(neisX,neisR,'-y')
    # plt.plot(waisX,waisR,'-b')


    # 绘图
    # plt.show()


    ######################################################################################################################


    # 二、下面是abaqus分析模块


    ######################################################################################################################


    # abaqus建模
    #引入abaqus模块建模
    # #导入abaqus模块
    # from abaqus import *
    # from abaqusConstants import *
    # from caeModules import *
    # import sketch

    Mdb()#创建新的模型数据库
    modelName = 'composite'
    m = mdb.Model(name=modelName)


    #在abaqus中绘制草图
    session.viewports['Viewport: 1'].view.fitView()
    s1 = m.ConstrainedSketch(name='composite-'+str('huojianfadongjiketi'), sheetSize=2*int(waisX[-1]))
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    constructionline=s1.ConstructionLine(point1=(waisX[0], 0), point2=(waisX[-1], 0))  # 根据外层X坐标自动调整轴线长度
    s1.assignCenterline(line=constructionline)
    s1.FixedConstraint(entity=constructionline)
    Vline=s1.ConstructionLine(point1=(0.0, -20.0), point2=(0.0, max(waisR)))
    s1.FixedConstraint(entity=Vline)


    ######################################################################################################################


    j=0
    for i in range(1,len(waisX)-2):
        if waisX[i+1-j]<tsX[0] and waisX[i-j]>min(Xnei[:]):
            if (waisX[i+1-j]-waisX[i-j])**2+(waisR[i+1-j]-waisR[i-j])**2<1:#距离太近的清除掉,防止后面abaqus的B样条畸形
                del(waisX[i-j])
                del(waisR[i-j])
                j += 1
        if waisX[i+1-j]>tsX[1] and waisX[i+1-j]<max(Xnei[:]):
            if (waisX[i+1-j]-waisX[i-j])**2+(waisR[i+1-j]-waisR[i-j])**2<1:
                del(waisX[i+1-j])
                del(waisR[i+1-j])
                j += 1


    ######################################################################################################################


    # 绘制轮廓草图
    Spoints1 = zip(neisX[:], neisR[:])#打包坐标点
    curve1 = s1.Spline(points=Spoints1)
    Spoints2 = zip(waisX[:], waisR[:])#打包坐标点
    curve2 = s1.Spline(points=Spoints2)
    line_cy5 = s1.Line(point1=(XddianZ[0], RddianZ[0]),
                       point2=(XddianZ[1], RddianZ[1]))  # 端点的存储没有区分纯螺旋与非纯螺旋
    line_cy6 = s1.Line(point1=(XddianY[0], RddianY[0]),
                       point2=(XddianY[1], RddianY[1]))

    # 生成三维实体模型
    session.viewports['Viewport: 1'].view.fitView()
    s1.setPrimaryObject(option=STANDALONE)
    s1.sketchOptions.setValues(constructionGeometry=ON)
    p = mdb.models['composite'].Part(name='Part-'+str('huojianfadongjiketi'), dimensionality=THREE_D,type=DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch=s1, angle=zhuanjiao, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()


    ######################################################################################################################


    # 导入装配模块
    a = mdb.models['composite'].rootAssembly
    #session.viewports['Viewport: 1'].setValues(displayedObject=a)#显示，#转换到装配体模块视图
    #session.viewports['Viewport: 1'].assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['composite'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)


    ######################################################################################################################


    # 调节显示精度
    p = mdb.models['composite'].parts['Part-'+str('huojianfadongjiketi')]
    p.setValues(geometryRefinement=EXTRA_FINE)
    a.Instance(name='Part-'+str('huojianfadongjiketi')+str(1), part=p, dependent=ON)


    ######################################################################################################################


    # 定义叠层上表面
    f,e=p.faces,p.edges
    lenn=len(waisX)/2
    lenn1=len(neisX)/2
    side1Faces22 = f.findAt(((waisX[lenn],waisR[lenn]*cos(0.5*zhuanjiao*todegree),waisR[lenn]*sin(0.5*zhuanjiao*todegree)), ))
    side1Faces22=side1Faces22+f.findAt(((waisX[lenn/2],waisR[lenn/2]*cos(0.5*zhuanjiao*todegree),waisR[lenn/2]*sin(0.5*zhuanjiao*todegree)), ))
    side1Faces22=side1Faces22+f.findAt(((waisX[4*lenn/3],waisR[4*lenn/3]*cos(0.5*zhuanjiao*todegree),waisR[4*lenn/3]*sin(0.5*zhuanjiao*todegree)), ))
    p.Surface(side1Faces=side1Faces22, name='Surf-top')#后面赋予属性可以用


    ######################################################################################################################


    # 切片疏密大小控制
    zuofengtou=[]
    youfengtou=[]
    # yuanshishuju=zip(X12,R12,alpha1,h2)#得到赋予角度的原始cadwind数据
    yuanshishuju=[[X12[i],R12[i],alpha1[i],h2[i]] for i in range(len(X12))]#得到赋予角度的原始cadwind数据
    pingjudistance=(max(X12)-abs(min(X12)))/len(X12)#得到这里的原始数据点平均距离
    # pingjudistance=1#得到这里的原始数据点平均距离
    # # 距离较近两点取舍
    shanchushu=[]
    for i in range(1,len(yuanshishuju)):
        distan=((yuanshishuju[i][0]-yuanshishuju[i-1][0])**2+(yuanshishuju[i][1]-yuanshishuju[i-1][1])**2)**0.5#两点的距离
        if distan<pingjudistance*0.3:
            shanchushu.append(yuanshishuju[i-1])
    # 距离较近两点取舍，删除距离较近的点
    for i in shanchushu[::-1]:
        yuanshishuju.remove(i)

    # 拆分左右封头
    for i in range(len(yuanshishuju)):
        if yuanshishuju[i][0]<=zuotongshen :
            # zuofengtou=yuanshishuju[0:i+1]
            zuofengtou.append(yuanshishuju[i])
            # break
    for i in range(len(yuanshishuju)):
    # for i in range(len(X12)):
        if yuanshishuju[i][0]>=youtongshen:
            # youfengtou=yuanshishuju[i:-1]
            youfengtou.append(yuanshishuju[i])
            # break


    ######################################################################################################################


    # 切分封头实体模型
    # 根据内表面的数据点，得到外表面的数据点
    def yitihuaqieFen(fengtou,huanxsl,luoxuansl):
        top_xvalue,top_yvalue,bot_xvalue,bot_yvalue=[],[],[],[]
        for i in range(len(fengtou)-1):
            k1=(fengtou[i+1][1]-fengtou[i][1])/(fengtou[i+1][0]-fengtou[i][0])
            if k1==0:
                top_xvalue.append(fengtou[i][0])
                top_yvalue.append((fengtou[i][1]+(3*luoxuansl)*fengtou[i][3]))
                bot_xvalue.append(fengtou[i][0])
                bot_yvalue.append((fengtou[i][1]-(0.2*luoxuansl)*fengtou[i][3]))
            else:
                k2=-1/k1
                if fengtou[i][0]<=zuotongshen and i<len(fengtou)-2:
                    theta1=arctan(k2)#弧度值
                    # print(theta1)
                    top_xvalue.append((fengtou[i+1][0]-(2*luoxuansl)*fengtou[i][3]*cos(theta1)))
                    top_yvalue.append((fengtou[i+1][1]-(2*luoxuansl)*fengtou[i][3]*sin(theta1)))
                    bot_xvalue.append((fengtou[i+1][0]+(0.2*luoxuansl)*fengtou[i][3]*cos(theta1)))
                    bot_yvalue.append((fengtou[i+1][1]+(0.2*luoxuansl)*fengtou[i][3]*sin(theta1)))
                elif fengtou[i][0]<=zuotongshen and i==len(fengtou)-2:
                    top_xvalue.append(zuotongshen)
                    top_yvalue.append((fengtou[i+1][1]+(huanxsl+2*luoxuansl)*fengtou[i][3]))
                    bot_xvalue.append(zuotongshen)
                    bot_yvalue.append((fengtou[i+1][1]-(0.2*luoxuansl)*fengtou[i][3]))
                elif fengtou[i][0] >= youtongshen:
                    theta1 = arctan(k2)  # 弧度值
                    if i==0:
                        top_xvalue.append(youtongshen)
                        top_yvalue.append((fengtou[i][1] + (huanxsl+2* luoxuansl) * fengtou[i][3] * sin(theta1)))
                        bot_xvalue.append(youtongshen)
                        bot_yvalue.append((fengtou[i][1] - (0.2 * luoxuansl) * fengtou[i][3] * sin(theta1)))
                    if i>0:
                        if theta1>=0:
                            top_xvalue.append(fengtou[i][0])
                            top_yvalue.append((fengtou[i][1] + (3* luoxuansl) * fengtou[i][3] * sin(theta1)))
                            bot_xvalue.append((fengtou[i][0]))
                            bot_yvalue.append((fengtou[i][1] - (0.2*luoxuansl) * fengtou[i][3] * sin(theta1)))
                        elif theta1<0:
                            top_xvalue.append(fengtou[i][0])
                            top_yvalue.append((fengtou[i][1] - (3 * luoxuansl) * fengtou[i][3] * sin(theta1)))
                            bot_xvalue.append((fengtou[i][0]))
                            bot_yvalue.append((fengtou[i][1] + (0.2*luoxuansl) * fengtou[i][3] * sin(theta1)))
        qiefenwailunkuo=zip(top_xvalue,top_yvalue)
        qiefenneilunkuo=zip(bot_xvalue,bot_yvalue)
        return qiefenwailunkuo,qiefenneilunkuo
    # 用厚度来分割
    zqiefenwailunkuo,zqiefenneilunkuo=yitihuaqieFen(zuofengtou,huanxsl,luoxuansl)
    yqiefenwailunkuo,yqiefenneilunkuo=yitihuaqieFen(youfengtou,huanxsl,luoxuansl)


    ######################################################################################################################


    # 进行切分侧面处理处理
    f, e, dpart= p.faces, p.edges, p.datums

    # 建立分割草图定位面
    pickFace=f.findAt(coordinates=((0.5*(youtongshen+zuotongshen),0.5*(max(neisR)+max(waisR)), 0.0)))#通过中点坐标得到曲面

    s = mdb.models['composite'].ConstrainedSketch(name='PartitionSketch',
        sheetSize=max(waisX)*3, gridSpacing=round(max(neisX)/10))
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    # 修建分割草图
    try:
        for k in range(0,len(zqiefenwailunkuo)):#建立分割草图线段
            line1=s.Line(point1=zqiefenwailunkuo[k],point2=zqiefenneilunkuo[k])
            # s.autoTrimCurve(curve1=line1, point1=zqiefenwailunkuo[k])#修建直线
            # line11=g.findAt(zqiefenneilunkuo[k])#寻找修建后的直线
            # s.autoTrimCurve(curve1=line11, point1=zqiefenneilunkuo[k])
        for k in range(0,len(yqiefenwailunkuo)):#建立分割草图线段
            line2=s.Line(point1=yqiefenwailunkuo[k],point2=yqiefenneilunkuo[k])
            # s.autoTrimCurve(curve1=line2, point1=yqiefenwailunkuo[k])
            # line22 = g.findAt(yqiefenneilunkuo[k])  # 寻找修建后的直线
            # s.autoTrimCurve(curve1=line22, point1=yqiefenneilunkuo[k])
    except:
        pass
        # print("AutoTrimCurve Error")
    # 定义一条边用于草图定位显示
    Axisline=p.DatumAxisByTwoPoint(point1=(waisX[0],0,0),point2=(waisX[-1],0,0))#创建轴线
    p.PartitionFaceBySketch(sketchUpEdge=dpart[Axisline.id], faces=pickFace,sketchOrientation=BOTTOM,sketch=s)#分割草图定位面


    ######################################################################################################################


    # 切分模型单元
    qiefenbian=[]
    quxianbian=[]

    f1, e1, d2 = p.faces, p.edges, p.datums
    for edge in e1:
        if abs(edge.pointOn[0][2])<1e-6 and edge.pointOn[0][1]>0 :#找到这里的z轴为0所在的面的所有的边，pointOn[0]元组好存储x,y,z,pointOn[0][2]得到z轴数据,有时候Z轴的坐标不是严格为0
            qiefenbian.append(edge)
            try:
                edge.getCurvature(parameter=0.5)#尝试有曲率的边，0.5为尝试数，只要不是0
                quxianbian.append(edge)
            except:#尝试曲率为0，则为直线，跳过
                pass

    quchu=[qiefenbian.remove(x) for x in quxianbian]#去除曲线边得到直线边

    refpoint=[]
    for i in range(len(qiefenbian)):
        dian=qiefenbian[i].pointOn#得到分割直线的点
        refpoint.append(dian[0])#得到参考点

    # 切分单元
    for point in refpoint:
        try:#尝试分割实体模型
            c=p.cells               
            pickedEdges =(e1.findAt(coordinates=point),)   
            p.PartitionCellBySweepEdge(sweepPath=e1.findAt(coordinates=(neisX[-1], neisR[-1]*cos(0.5*zhuanjiao*todegree), neisR[-1]*sin(0.5*zhuanjiao*todegree))), cells=c, edges=pickedEdges)        
        except:#尝试切分模型，遇到错误信息
            pass


    ######################################################################################################################


    # 为了在筒身中部添加约束，保证两边极孔的受力形式一样，需将筒身单元切分为两半
    pointts=(youtongshen/2, 0.0, 0.0)
    Planepoint=p.DatumPointByCoordinate(coords=pointts)#创建过平面的点
    Axisline=p.DatumAxisByTwoPoint(point1=(waisX[0],0,0),  point2=(waisX[-1],0,0))#创建轴线
    d_part = p.datums
    PlaneByPoint=p.DatumPlaneByPointNormal(point=d_part[Planepoint.id], normal=d_part[Axisline.id])
    # 分割中间cell
    pickedCells = c.findAt((youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
    , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180), ))#得到筒身cell
    p.PartitionCellByDatumPlane(datumPlane=d_part[PlaneByPoint.id], cells=pickedCells)


    ######################################################################################################################


    # 得到各个小片区的中点坐标，为后面建立cell赋属性用
    # 初始化相关的曲面数据，后面需要
    Masterface=[]
    # loadsurfacecenterpoint=[]
    f_all=p.faces
    face_coordinate=[]
    face_coordinate1=[]
    face_coordinate2=[]
    face_coordinate3=[]

    for face in f_all:
        centercoordinate = face.getCentroid()  # 将曲面一个一个的输入判断
        # print (centercoordinate)
        if centercoordinate[0][2] == 0.0 and centercoordinate[0][1] > 0.0:  # 根据面的中心点，得到z为0的面
            Masterface.append(centercoordinate)
            if centercoordinate[0][0] < tsX[0]:
                face_coordinate1.append(centercoordinate[0])
            if centercoordinate[0][0] >= tsX[0] and centercoordinate[0][0] <= tsX[1]:
                face_coordinate2.append(centercoordinate[0])
            if centercoordinate[0][0] > tsX[1]:
                face_coordinate3.append(centercoordinate[0])
    # face_coordinate1.sort(key=(lambda x: x[0]))#这里的角度排序有问题，因为X不是单调递增的，改动如下
    face_coordinate1.sort(key=(lambda x: x[1]))  # 这里的角度排序有问题，因为X不是单调递增的，封头按半径排序
    face_coordinate2.sort(key=(lambda x: x[0]))  # 这里的角度排序有问题，因为X不是单调递增的，筒身按X排序
    face_coordinate3.sort(key=(lambda x: -x[1]))  # 这里的角度排序有问题，因为X不是单调递增的，封头按半径排序
    # print(len(face_coordinate1))
    # print(face_coordinate1)
    # print(len(face_coordinate2))
    # print(face_coordinate2)
    # print(len(face_coordinate3))
    # print(face_coordinate3)
    for i in range(len(face_coordinate1) + len(face_coordinate2) + len(face_coordinate3)):
        if i < len(face_coordinate1):
            face_coordinate.append(face_coordinate1[i])
        if i < (len(face_coordinate1) + len(face_coordinate2)) and i >= len(face_coordinate1):
            face_coordinate.append(face_coordinate2[i - ((len(face_coordinate1) + len(face_coordinate2)))])
        if i >= (len(face_coordinate1) + len(face_coordinate2)) and i < (
                len(face_coordinate1) + len(face_coordinate2) + len(face_coordinate3)):
            face_coordinate.append(face_coordinate3[i - ((len(face_coordinate1) + len(face_coordinate2)))])
    # print(face_coordinate)


    ######################################################################################################################


    #创建复合材料工程常数
    mat_comp=m.Material(name='Material-composite')
    mat_comp.Elastic(type=ENGINEERING_CONSTANTS, 
        table=(engineeringdata, ))
    m.HomogeneousSolidSection(name='Section-composite', material='Material-composite', thickness=None)

    list_angle=[]
    for i in range(len(zuofengtou)):#得到cadwind过来的缠绕角数据
        list_angle.append(zuofengtou[i][2])
    list_angle.append(0)
    for i in range(len(youfengtou)):
        list_angle.append(youfengtou[i][2])

    count=0
    tongshendian=[]#分别建立筒身和封头的set集
    for k in face_coordinate:#从所有的z轴为0的面的中点内得到K数据
        c=p.cells      
        count=count+1
        pickcells=c.findAt(((k[0], k[1], k[2]),))
        p.Set(name='Set'+'-'+str(count),cells=pickcells)#创建集合
        if count==len(zuofengtou)+1 or  count==len(zuofengtou)+2:
            tongshendian.append((k[0], k[1], k[2]),)#创建筒身数据点


    ######################################################################################################################


    # 封头赋予铺层角度属性
    for k in range(1,count+1):
        if k!=len(face_coordinate1)+1 or k!=len(face_coordinate1)+2:#左封头铺层
            layupOrientation = None
            region1 = p.sets['Set-'+str(k)]
            normalAxisRegion = p.surfaces['Surf-top']
            compositeLayup = p.CompositeLayup(
                name='CompositeLayup-'+str(k), description='', elementType=SOLID,
                symmetric=False)
            # compositeLayup = p.CompositeLayup(
            #     name='CompositeLayup-'+str(k), description='', elementType=CONTINUUM_SHELL,
            #     symmetric=False)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
                poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
                useDensity=OFF)
            for i in range(1,luoxuansl+1):
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*i-1), region=region1, 
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*i), region=region1, 
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=1.0, 
                    orientationType=SPECIFY_ORIENT, orientationValue=-1*list_angle[k-1], 
                    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
            compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None, 
                additionalRotationType=ROTATION_NONE, angle=0.0, 
                additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3, 
                normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion, 
                normalAxisDirection=AXIS_3, flipNormalDirection=False, 
                primaryAxisDefinition=VECTOR, primaryAxisVector=(1, 0.0, 0.0), 
               primaryAxisDirection=AXIS_1, flipPrimaryDirection=False) 
        else:
            pass


    ######################################################################################################################


    ###筒身赋予角度
    layupOrientation = None
    for i in range(1,3):
        region1 = p.sets['Set-'+str(len(face_coordinate1)+i)]
        normalAxisRegion = p.surfaces['Surf-top']
        compositeLayup = p.CompositeLayup(
            name='CompositeLayup-'+str(len(face_coordinate1)+i), description='', elementType=SOLID,
            symmetric=False)
        # compositeLayup = p.CompositeLayup(
        #     name='CompositeLayup-'+str(len(face_coordinate1)+i), description='', elementType=CONTINUUM_SHELL,
        #     symmetric=False)
        compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON,
            poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT,
            useDensity=OFF)
        ply_count=0
        for x in puceng:
            if x=='HelixLoop':
                ply_count=ply_count+1
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*ply_count-1), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=lxdt,
                    orientationType=SPECIFY_ORIENT, orientationValue=luoxuanjiao,
                    additionalRotationType=ROTATION_NONE, additionalRotationField='',
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-' + str(2*ply_count), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=lxdt,
                    orientationType=SPECIFY_ORIENT, orientationValue=-luoxuanjiao,
                    additionalRotationType=ROTATION_NONE, additionalRotationField='',
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
            elif x=='Hoop':
                ply_count=ply_count+1
                compositeLayup.CompositePly(suppressed=False, plyName='Ply-'+str(2*ply_count-1), region=region1,
                    material='Material-composite', thicknessType=SPECIFY_THICKNESS, thickness=hxdt,
                    orientationType=SPECIFY_ORIENT, orientationValue=90,
                    additionalRotationType=ROTATION_NONE, additionalRotationField='',
                    axis=AXIS_3, angle=0.0, numIntPoints=3)
        compositeLayup.ReferenceOrientation(orientationType=DISCRETE, localCsys=None,
            additionalRotationType=ROTATION_NONE, angle=0.0,
            additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3,
            normalAxisDefinition=SURFACE, normalAxisRegion=normalAxisRegion,
            normalAxisDirection=AXIS_3, flipNormalDirection=False,
            primaryAxisDefinition=VECTOR, primaryAxisVector=(1, 0.0, 0.0),
           primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)


    ######################################################################################################################


    # 创建step
    mdb.models['composite'].StaticStep(name='Step-1', previous='Initial',
        maxNumInc=300, initialInc=0.05, maxInc=0.1)
    mdb.models['composite'].steps['Step-1'].setValues(nlgeom=ON)

    ############################################################
    # 创建场变量输出
    for i in range(2,len(face_coordinate)+2):
        mdb.models['composite'].FieldOutputRequest(name='F-Output-'+str(i),
            createStepName='Step-1', variables=('S', 'E', 'LE'), layupNames=(
            'Part-huojianfadongjiketi1.CompositeLayup-'+str(i-1), ),
            layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True,
            outputAtPlyBottom=False, rebar=EXCLUDE)


    ######################################################################################################################


    # 通过旋转侧面的中点坐标创建柱坐标
    # DatumcsysCylinder=a.DatumCsysByThreePoints(\
    #     point1=(face_coordinate2[0][0],face_coordinate2[0][1]*math.cos(zhuanjiao/2*math.pi/180),\
    #                                 face_coordinate2[0][1]*math.sin(zhuanjiao/2*math.pi/180)), name='Datum csys-Cylinder',
    #     coordSysType=CYLINDRICAL, origin=(face_coordinate2[0][0], 0.0, 0.0),
    #     point2=face_coordinate2[0])

    DatumcsysCylinder=a.DatumCsysByThreePoints(\
        point1=(youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
        , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180)), name='Datum csys-Cylinder1',
        coordSysType=CYLINDRICAL, origin=(youtongshen/2, 0.0, 0.0),
        point2=(youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
        , 0))


    ######################################################################################################################


    # 创建参考点
    rp=a.ReferencePoint(point=pointts)
    ReferencePoints=a.referencePoints#assembly层次的参考点存储仓库
    region1=a.Set(referencePoints=[ReferencePoints[rp.id],], name='m_Set-1')
    f1 = a.instances['Part-huojianfadongjiketi1'].faces
    faceszhong = f1.findAt(((youtongshen/2,0.5*(max(neisR)+max(waisR))*math.cos(zhuanjiao/2*math.pi/180)\
                                                                            , 0.5*(max(neisR)+max(waisR))*math.sin(zhuanjiao/2*math.pi/180)), ))
    # region2=regionToolset.Region(faces=faceszhong)
    # datum =a.datums[DatumcsysCylinder.id]#从坐标系仓库中获取到柱坐标
    # # 建立耦合约束
    # mdb.models['composite'].Coupling(name='Constraint-1', controlPoint=region1,
    #     surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
    #     localCsys=datum, u1=OFF, u2=OFF, u3=ON, ur1=ON, ur2=ON, ur3=ON)


    ######################################################################################################################


    # 找到主面的中点坐标
    Slaveface=[]
    for i in range(len(Masterface)):
        abc1=Masterface[i][0][0]
        abc2=Masterface[i][0][1]*math.cos(math.pi*zhuanjiao/180)
        abc3=Masterface[i][0][1]*math.sin(math.pi*zhuanjiao/180)
        Slaveface.append(((abc1,abc2,abc3), ))
    # 创建主从面
    s1=a.instances['Part-huojianfadongjiketi1'].faces
    side1Faces1=s1.findAt(Masterface[0])
    for i in range(1,len(Masterface)):
        faces1=s1.findAt(Masterface[i])
        side1Faces1=side1Faces1+faces1#遍历读取面
    regionMaterface=a.Surface(side1Faces=side1Faces1, name='m_Surf-1')
    side1Faces2=s1.findAt(Slaveface[0])
    for i in range(1,len(Slaveface)):
        faces2 = s1.findAt(Slaveface[i])
        side1Faces2+=faces2#遍历读取面
    regionSlaveface=a.Surface(side2Faces=side1Faces2, name='s_Surf-1')
    # 创建参考带点，循环对称用
    r1 = a.ReferencePoint(point=(0.0, 0.0, 0.0))
    r1points=a.referencePoints
    regionreferencepoint1=regionToolset.Region(referencePoints=[r1points[r1.id]])#该种做法参见《abaquspython二次开发攻略》203页
    # region3=regionToolset.Region(referencePoints=refPoints1)
    r1 = a.ReferencePoint(point=(100.0, 0.0, 0.0))
    rpoints=a.referencePoints
    # refPoints1=(r1[2], )
    regionreferencepoint2=regionToolset.Region(referencePoints=[r1points[r1.id]])#该种做法参见《abaquspython二次开发攻略》203页
    # region4=regionToolset.Region(referencePoints=refPoints1)
    # 建立循环对称相互作用
    mdb.models['composite'].CyclicSymmetry(name='Int-1', createStepName='Initial',
        master=regionMaterface, slave=regionSlaveface, axisPoint1=regionreferencepoint1, axisPoint2=regionreferencepoint2,
        positionToleranceMethod=COMPUTED_TOLERANCE, positionTolerance=0.0,
        adjustTie=True, repetitiveSectors=360/int(zhuanjiao),
        extractedNodalDiameter=ALL_NODAL_DIAMETER, excitationNodalDiameter=0)
    #: The interaction "Int-1" has been created.


    ######################################################################################################################


    # 2020.10.23
    # 由前面的内轮廓数据点坐标，得到内轮廓曲面旋转后的每个小面的中点坐标
    # 将建立内表面的
    # 2020.9.26改，给每个cell找到面，筛选出内外表面，进而得到内表面
    # 1.遍历每一个cell
        ##1.1得到cell的每1个面
        ##1.2根据曲面的曲率判定上下弧面,内外弧面含有两个曲率半径，其他要么是平面要么其他曲面
    #2.将得到的内外表面上数据点按径向距离区分内外表面
    # 建立内表面的几何surfaces
    for iijj  in range((len(face_coordinate))):#cell的数量是和侧面中点数是一样的
        cellface=p.cells[iijj].getFaces()#得到当前cell的6个表面
        facex = []
        neiwai=[]
        for ih in cellface:
            facepoint0 = p.faces[ih].pointOn[0]  # 得到当前cell当前编号face的随机面内坐标,不能用中点，因为曲面的中点不一定在面上
            try:
                p.faces[ih].getCurvature(0, 0)["curvature1"] != 0 or p.faces[ih].getCurvature(0, 0)["curvature2"] != 0
                qulv1 = p.faces[ih].getCurvature(0, 0)["curvature1"]
                # print(qulv1)
                qulv2 = p.faces[ih].getCurvature(0, 0)["curvature2"]
                # print(qulv2)
                if qulv1  != 0 and qulv2 != 0:
                    neiwai.append((facepoint0,ih))
                    facex.append(p.faces[ih].pointOn[0][0])
                    # print(neiwai)
            except:
                pass
        if len(facex)==2:
            print("waibiaomian")
            if (facex[0] + facex[1]) / 2 >= zuotongshen and (facex[0] + facex[1]) / 2 <= youtongshen:  # 先判断是不是筒身，筒身的中点连判断，因为筒身这里的判断有嗲小问题不用这个的话
                print(facex)
                neiwai1=[neiwai[i][0] for i in range(2)]
                neiwai1.sort(key=(lambda x: (x[1] ** 2 + x[2] ** 2)))  # 按旋转半径排序
                # print(neiwai)
                # print(neiwai1)
                # neisurfaces.append((neiwai1[0][0],))
                # neiwaiSurfaces.append(neiwai)  # 存储内外表面
                break #得到即退出
    StackDirectionSurface = neiwai1[-1]
    neiwai.sort(key=(lambda x: (x[0][1] ** 2 + x[0][2] ** 2)))  # 按旋转半径排序
    neibiaomian = p.faces[neiwai[0][1]]
    neibiaomian1 = neibiaomian.getFacesByFaceAngle(10)

    #根据需求截断x
    neisurfaces=[]
    for i in range(len(neibiaomian1)):
        if neibiaomian1[i].pointOn[0][0]<=tschangdu+jinshunei:
            neisurfaces.append(neibiaomian1[i].pointOn)

    side1Faces3= s1.findAt(neisurfaces[0])#通过中面寻找到曲面
    for i in range(1,len(neisurfaces)):
        faces3=s1.findAt(neisurfaces[i])
        side1Faces3=side1Faces3+faces3#遍历读取面

    region5=a.Surface(side1Faces=side1Faces3, name='Surf-bottom')
    mdb.models['composite'].Pressure(name='Load-1', createStepName='Step-1',
        region=region5, distributionType=UNIFORM, field='', magnitude=pressure,
        amplitude=UNSET)


    ######################################################################################################################


    # 建立边界条件1,对旋转周期对称约束进行边界约束
    regionBC = regionToolset.Region(faces=side1Faces1)
    datum = a.datums[DatumcsysCylinder.id]
    mdb.models['composite'].YsymmBC(name='BC-1', createStepName='Initial',
        region=regionBC, localCsys=datum)

    # 边界条件2，对筒身中部这里的的轴向进行限定
    regiontshen= regionToolset.Region(faces=faceszhong)
    mdb.models['composite'].DisplacementBC(name='BC-2', createStepName='Initial',
        region=regiontshen, u1=UNSET, u2=UNSET, u3=SET, ur1=SET, ur2=SET, ur3=SET,
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)


    ######################################################################################################################


    # 2020.10.22
    # 种子点分布控制，厚度方向一个网格
    p.deleteMesh()
    # 思路：
    # 1.由之前旋转侧面的中点数量坐标得到cell数量
    # 2.得到每一个cell的边
    # 3.找到cell的各边单独控制网格数量
    # 建立cell的边数据
    cellzuoyouzong=[]#所有单元左右边
    cellshangxiazong=[]#所有单元左右边
    cellcurveszong=[]#所有单元左右边
    for iijj  in range((len(face_coordinate))):#cell的数量是和侧面中点数是一样的
        celledges=p.cells[iijj].getEdges()#得到当前cell的12个边
        edgepoint0=[]
        edgepoint1=[]
        hubian=[]
        pp=0
        for ih in celledges:
            if p.edges[ih].pointOn[0][2]<1e-6:
                edgepoint0.append((p.edges[ih].pointOn[0],ih))#得到当前cell当前编号edge的线上随机坐标
                # print('edgepoint0')
            if (zhuanjiao-1e-6)<math.atan(p.edges[ih].pointOn[0][2]/p.edges[ih].pointOn[0][1])*180/math.pi<=(zhuanjiao+1e-6):
                edgepoint1.append((p.edges[ih].pointOn[0],ih))#得到当前cell当前编号edge的线上随机坐标
                # print('edgepoint1')
            if 0<math.atan(p.edges[ih].pointOn[0][2]/p.edges[ih].pointOn[0][1])/math.pi*180<zhuanjiao:
                hubian.append((p.edges[ih].pointOn[0],ih))
    # 2020.11.1
        if len(edgepoint0)>0:
            Xx0=sum([edgepoint0[i][0][0] for i in range(len(edgepoint0))])/len(edgepoint0)
            aa=edgepoint0
            # print(edgepoint0)
            # print(Xx0)
        if len(edgepoint1)>0:
            Xx1=sum([edgepoint1[i][0][0] for i in range(len(edgepoint1))])/len(edgepoint1)
            bb=edgepoint1
        Xx0=(Xx0+Xx1)/2
        if tsX[0] < Xx0 < tsX[1]:  # 筒身
            if len(aa)>0:
                dd=aa
                dd.sort(key=(lambda x: (x[0][0] )))  # 旋转半径排序
                point = p.edges[dd[len(dd)/2][1]]
                edges = point.getEdgesByEdgeAngle(20)
                length = []
                j = 0
                if j == 0:
                    for i in range(len(edges)):
                        length.append(edges[i].getSize(printResults=False))
                    cc = sorted(length)[0:-len(length) / 2 + 1]
                    minsizes = min(cc)
                    Avlength = sum(cc) / len(edges)
                    j += 1
                # print(1)
                # print(Avlength)
                # aa.sort(key=(lambda x: (x[0][1] ** 2 + x[0][2] ** 2)))  # 旋转半径排序
                aa.sort(key=(lambda x: (x[0][0])))  # 旋转半径排序
                # print(aa)
                partedges = p.edges
                for i in range(len(aa)):
                    # 上下边
                    # if i <=0 or i >=3:#内表面边，上下边
                    if i ==0 or i ==2:#内表面边，上下边
                        EdgeLength = p.edges[aa[i][1]].getSize(printResults=False)  # 获取边长度
                        if EdgeLength<Avlength:
                            number=1
                        elif EdgeLength>Avlength:
                            number = int((EdgeLength/Avlength)*5)+1
                        upANDdownEdges = partedges.findAt((aa[i][0],))
                        p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 上下边布种子
                    # 左右边
                    else:
                        cellzuoyouzong.append((aa[i][0],))  # 直边即是左右
                        leftANDraghtEdges = partedges.findAt((aa[i][0],))
                        p.seedEdgeByNumber(edges=leftANDraghtEdges, number=1)  # 左右边布种子
            if len(bb)>0:
                # bb.sort(key=(lambda x: (x[0][1] ** 2 + x[0][2] ** 2)))  # 旋转半径排序
                bb.sort(key=(lambda x: (x[0][0])))  # 旋转半径排序
                # print(bb)
                partedges = p.edges
                for i in range(len(bb)):
                    # 上下边
                    # if i <=0 or i >=3:#内表面边，上下边
                    if i ==1 or i ==2:#内表面边，上下边
                        EdgeLength = p.edges[bb[i][1]].getSize(printResults=False)  # 获取边长度
                        if EdgeLength < Avlength:
                            number = 1
                        elif EdgeLength > Avlength:
                            number = int((EdgeLength/Avlength)*10) + 1
                        upANDdownEdges = partedges.findAt((bb[i][0],))
                        p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 上下边布种子
                    # 左右边
                    else:
                        cellzuoyouzong.append((bb[i][0],))  # 直边即是左右
                        leftANDraghtEdges = partedges.findAt((bb[i][0],))
                        p.seedEdgeByNumber(edges=leftANDraghtEdges, number=1)  # 左右边布种子
        elif tsX[0] > Xx0 or Xx0> tsX[1]:  # 封头
            if len(aa)>0:
                for i in range(len(aa)):
                    partedges = p.edges
                    try:  # 先尝试曲边寻找
                        p.edges[aa[i][1]].getCurvature(0.5)
                        cellshangxiazong.append((aa[i][0],))  # 曲边即是上下
                        EdgeLength = p.edges[bb[i][1]].getSize(printResults=False)  # 获取边长度
                        if EdgeLength < Avlength:
                            number = 1
                        elif EdgeLength > Avlength:
                            number = int(EdgeLength / Avlength*3) + 1
                        upANDdownEdges = partedges.findAt((aa[i][0],))
                        p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 上下边布种子
                    except:
                        cellzuoyouzong.append((aa[i][0],))  # 直边即是左右
                        leftANDrightEdges = partedges.findAt((aa[i][0],))
                        p.seedEdgeByNumber(edges=leftANDrightEdges, number=1)  # 左右边布种子
            if len(bb)>0:
                for i in range(len(bb)):
                    partedges = p.edges
                    try:  # 先尝试曲边寻找
                        p.edges[bb[i][1]].getCurvature(0.5)
                        cellshangxiazong.append((bb[i][0],))  # 曲边即是上下
                        EdgeLength = p.edges[bb[i][1]].getSize(printResults=False)  # 获取边长度
                        if EdgeLength < Avlength:
                            number = 1
                        elif EdgeLength > Avlength:
                            number = int(EdgeLength / Avlength*3) + 1
                        upANDdownEdges = partedges.findAt((bb[i][0],))
                        p.seedEdgeByNumber(edges=upANDdownEdges, number=number)  # 上下边布种子
                    except:
                        cellzuoyouzong.append((bb[i][0],))  # 直边即是左右
                        leftANDrightEdges = partedges.findAt((bb[i][0],))
                        p.seedEdgeByNumber(edges=leftANDrightEdges, number=1)  # 左右边布种子
        # 弧边
        i=0
        if i ==0:
            for i in range(len(hubian)):
                cellcurveszong.append((hubian[i][0],))
                number = int(p.edges[hubian[i][1]].getSize(printResults=False) / 5)  # 获取边长度并除以3得到单元布置数量初值
                if number == 0:
                    number = 1
                CurvesEdges = partedges.findAt((hubian[i][0],))
                p.seedEdgeByNumber(edges=CurvesEdges, number=number)  # 左右边布种子
            i+=1
    # 网格绘制
    p.generateMesh()


    # 为各边不同的网格控制方式生成相应的set集
    # 厚度方向网格数量控制
    partedges = p.edges
    leftANDraghtEdges = partedges.findAt(cellzuoyouzong[0])
    for ii in range(1,len(cellzuoyouzong)):
        Edge= partedges.findAt(cellzuoyouzong[ii])
        leftANDraghtEdges =leftANDraghtEdges+Edge
    # 创建左右边集合
    p.Set(name='Set'+'-'+'HouduEdges',edges=leftANDraghtEdges)
    # p.seedEdgeByNumber(edges=leftANDraghtEdges,number=1,constraint=FINER)
    # 经线方向上网格数量控制
    upANDdownEdges = partedges.findAt(cellshangxiazong[0])
    for ii in range(1,len(cellshangxiazong)):
        Edge= partedges.findAt(cellshangxiazong[ii])
        upANDdownEdges =upANDdownEdges+Edge
    # 创建上下边集合
    p.Set(name='Set'+'-'+'JingxianEdges',edges=upANDdownEdges)
    # p.seedEdgeByNumber(edges=upANDdownEdges,number=5,constraint=FINER)#非筒身布种子
    # 旋转弧线上的网格数量控制
    CurvesEdges = partedges.findAt(cellcurveszong[0])
    for ii in range(1,len(cellcurveszong)):
        Edge= partedges.findAt(cellcurveszong[ii])
        CurvesEdges =CurvesEdges+Edge
    # 创建弧线边集合
    p.Set(name='Set'+'-'+'HuxianEdges',edges=CurvesEdges)

    # 单元类型选择
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD,
        kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
        hourglassControl=DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)

    # elemType1 = mesh.ElemType(elemCode=SC8R, elemLibrary=STANDARD,
    #     secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
    # elemType2 = mesh.ElemType(elemCode=SC6R, elemLibrary=STANDARD)
    # elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)



    # 选择cell
    c = p.cells
    Regionszhiding=c.findAt((face_coordinate[0],))
    for k in range(1,len(face_coordinate)):#因为前面只用到了z=0的侧面，所以可以通过该面的中点得到cells
        Regionszhiding1=c.findAt((face_coordinate[k],))
        Regionszhiding=Regionszhiding+Regionszhiding1
    cell=p.Set(name='Set' + '-' + 'cell', cells=Regionszhiding)
    p.setElementType(regions=cell, elemTypes=(elemType1, elemType2,elemType3))

    # 转换叠层上下表面
    f = p.faces
    p.assignStackDirection(referenceRegion=f.findAt(coordinates=StackDirectionSurface), cells=Regionszhiding)


    ######################################################################################################################



    # 切换到JOB模块
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    session.viewports['Viewport: 1'].view.fitView()

    # 创建一个job
    mdb.Job(name='Job-1-composite-huojianfadongjiketi', model='composite',
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)

    # 提交job
    # mdb.jobs['Job-1-composite-huojianfadongjiketi'].submit(consistencyChecking=OFF)


    ######################################################################################################################








