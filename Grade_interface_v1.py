####execfile('C:/Scratch/Finger/Grade_interface_v1.py')

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import*
from regionToolset import *
import regionToolset
import math
import mesh
import random
import numpy as np
import os
from material import createMaterialFromDataString
import sys
sys.path.insert(15, r'c:/Programs/SIMULIA/CAE/2017/win_b64/code/python2.7/lib/abaqus_plugins/stlExport')
#import stlExport_kernel

execfile('C:\\Scratch\\Finger\\AbaqusScriptFunc.py')
os.chdir(r"C:\\Scratch\\Finger")
subfolder = '\\Finger'
picfolder = 'C:\\Scratch\\Finger\\pictures'
secondfolder = 'C:\\Scratch' + subfolder

type_lab = ['Finger_v1_grey']
strainperc = - 10.

divs = 1. ## THIS NUMBER WILL SUBDIVIDE THE MESH TO OBTAIN MORE ELEMENTS WITHOUT CHANGING ANY DIMENSION

#DONT TOUCH THESE:
perc_vec = 1.0 #no label means it was -33p
resvecval =1.0
tol=0.0001
CPUSS = 6


# MATERIAL PROPERTY EQUATIONS
crackBC = 0
punish = 1.0
D_fin = 0.05

pointsnum = 500

alpha_soft_exp =    3.00629914
mu_soft_exp =    0.266108723
#alpha_soft_exp =    5.70951890
#mu_soft_exp =    0.331586536
D1 = 2/  ( 2*mu_soft_exp*(1+0.485)/(3*(1-2*0.485)) )  
eps_pl_soft_HE =     [0,0.1]
sig_pl_soft_HE =    [2.5160,3.0]

############
Ym_soft_LIU = 1.0

####### HARD MATERIAL AND SOFT MATERIAL, LIU CONSTANTS
fit_eps_ult = np.array(([0.697891058966723,9.21806674451016,0.452419504617438]))

#fit_A = np.array(([52.3795986952229,1.99250654729608,0.752050404057717]))
order_A = 2
theta_A = np.transpose(np.array(([[0.456238775513901,23.6763999666868,-16.7323196046030]])))

fit_A2 = np.array(([77.0329055702472,1,-23.6431115795406]))

fit_A20_1 = np.array(([66.3995423700196,1.69995902315442,0.694969706195360]))
fit_A40_1 = np.array(([14.9518547422338,0.882217235476243,0.694969706195360]))
fit_A20_2 = np.array(([51.5131099238831,2.12443880688038,3.31313960647499]))
fit_A40_2 = np.array(([67.9800142205931,1.27242355224625,-14.0151709781162]))

fit_alpha1 = np.array(([8.44207355238249,0.726412345470195,2.03536795433075]))
fit_alpha2 = np.array(([116.521386242516,1,-18.6464503889803]))
fit_alpha2 = np.array(([257.347350822932,0.2886254285243492,-157.067927248874]))

fit_alpha20_1 = np.array(([8.44207355238249,0.726412345470195,2.03536795433075]))
fit_alpha40_1 = np.array(([642.599923901834,3.36073599046238,2.03536795433075]))
fit_alpha20_2 = np.array(([226.986227089168,0.332590530310889,-128.243311766979]))
fit_alpha40_2 = np.array(([4778.62401712850,0.0157265230658667,-4678.65892877455]))

fit_ratio = np.array(([0.892459969593177,7.52046968494233,1.00371385070869]))

#############################################

cont  = 0

topcol = [0,102,162]
botcol = [195,49,47]



fflipp = 0
#Mdb()
cont  = 0
for II in range(0, len(type_lab)):
    modelo =  type_lab[II] 
    ## Time to begin
    print ' ' 
    print '        ' + modelo
    print 'DO YOU WRESTLE WITH DREAMS?'
    
    filename = secondfolder+ '\\abbas\\' + modelo + '.abba'
    filename_hulshult = secondfolder+ '\\abbas\\' + type_lab[0] + '.hulshult'
    ff = open(filename_hulshult,'r')
    lines = ff.readlines()
    ff.close()
    #modelo =  type_lab[II] + 'D'
    
    voxY = float(lines[0])
    voxX = float(lines[1])
    voxZ = float(lines[2])
    
    H = int(lines[3])
    W = int(lines[4])
    L = int(lines[5])
    
    Hr = H
    Wr = W
    Lr = L
        
    ff = open(filename,'r')
    lines = ff.readlines()
    ff.close()
    
    lin_num = len(lines)
    pmat = np.zeros((lin_num,12))
    CCC = np.array(())
    for ii in range(0,lin_num):
        aux = [float(jj) for jj in lines[ii].replace('(',' ').replace(')',' ').replace(',',' ')[0:len(lines[ii].replace('(',' ').replace(')',' ').replace(',',' '))].split()]
        CCC = np.concatenate((CCC, aux))
    
    CCC = CCC.reshape((H,W,L),order='F')[::-1,:,:]
    CCC_KILL = np.squeeze((CCC == 0).astype(int))
    
    CCC = np.repeat(np.repeat(CCC,divs,axis=0),divs,axis=1)# CASE SPECIFIC!!!!!!!!
    CCC = np.repeat(CCC,divs,axis=2)# CASE SPECIFIC!!!!!!!!
    CCC_KILL = np.repeat(np.repeat(CCC_KILL,divs,axis=0),divs,axis=1)# CASE SPECIFIC!!!!!!!!
    CCC_KILL = np.repeat(CCC_KILL,divs,axis=2)# CASE SPECIFIC!!!!!!!!
    voxX=voxX/divs# CASE SPECIFIC!!!!!!!!
    voxY=voxY/divs# CASE SPECIFIC!!!!!!!!
    voxZ=voxZ/divs# CASE SPECIFIC!!!!!!!!
    Hr=Hr*divs# CASE SPECIFIC!!!!!!!!
    Wr=Wr*divs# CASE SPECIFIC!!!!!!!!
    Lr=Lr*divs# CASE SPECIFIC!!!!!!!!
    
    if fflipp == 1:
        CCC = np.concatenate((CCC[::-1,:,:],CCC),axis=0)
        CCC_KILL = np.concatenate((CCC_KILL[::-1,:,:],CCC_KILL),axis=0)
        Hr = Hr*2
        voxZ = voxZ*2.
        #Wr = Wr*2
        #Lr = Lr*2
    
    tot_props= np.unique(CCC).size 
    prop_list = np.unique(CCC)
    
    
    propMatC=np.zeros(tot_props)
    
    for ii in range(0,tot_props):
        propMatC[ii] = np.sum(np.squeeze((CCC == prop_list[ii]).astype(int)))
    
    el_lab_size = int(propMatC.max())
    materialMat= np.zeros((tot_props,el_lab_size))
    
    rango=range(0,int(Hr))
    rangorevY=rango[::-1]
    rango=range(0,int(Wr))
    rangorevX=rango[::-1]
    rango=range(0,int(Lr))
    rangorevZ=rango[::-1]
    #oopsi
    ccount = 1
    propMatC=np.zeros(tot_props)
    lab23d_H = np.zeros(int(Wr)*int(Hr)*int(Lr)+1) #<-- give the label of an element, it will give you its indeces
    lab23d_W = np.zeros(int(Wr)*int(Hr)*int(Lr)+1)
    lab23d_L = np.zeros(int(Wr)*int(Hr)*int(Lr)+1)
    for ii in range(0,int(Hr)):
        for  kk in rangorevZ:
            for  jj in rangorevX:
                lab23d_H[ccount] = ii
                lab23d_W[ccount] = jj
                lab23d_L[ccount] = kk
                for mm in range(0,int(tot_props)): #mm marks 0 if is soft and 1 if its hard, so we need to compare with CCC_B and CCC_A respectively
                    if CCC[ii,jj,kk]==prop_list[mm]:
                        materialMat[mm,int(propMatC[mm])]=ccount
                        propMatC[mm]=propMatC[mm]+1
                        ccount = ccount + 1
    
    materialMat= materialMat.astype(int)
    print 'DO YOU CONTEND WITH SHADOWS?'
    
    #if it == 0:
    
    mdb.Model(name=modelo, modelType=STANDARD_EXPLICIT)
    mdb.models[modelo].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    P1=    (0.0,      0.0)
    P2=    (Wr*voxX,    0.0)
    P3=    (Wr*voxX,    Hr*voxY)
    P4=    (0.0,      Hr*voxY)
    mdb.models[modelo].sketches['__profile__'].Line(point1=P1, point2=P2)
    mdb.models[modelo].sketches['__profile__'].Line(point1=P2, point2=P3)
    mdb.models[modelo].sketches['__profile__'].Line(point1=P3, point2=P4)
    mdb.models[modelo].sketches['__profile__'].Line(point1=P4, point2=P1)
    mdb.models[modelo].Part(dimensionality=THREE_D, name='block', type=DEFORMABLE_BODY)
    mdb.models[modelo].parts['block'].BaseSolidExtrude(sketch=mdb.models[modelo].sketches['__profile__'],depth=Lr*voxZ)
    mdb.models[modelo].parts['block'].seedEdgeBySize(
        edges=mdb.models[modelo].parts['block'].edges.getByBoundingBox(
        xMin=-tol, yMin=-tol, zMin=-tol, xMax=Wr*voxX+tol, yMax=tol, zMax=tol) , 
        size=voxX, constraint=FIXED)
    mdb.models[modelo].parts['block'].seedEdgeBySize(
        edges=mdb.models[modelo].parts['block'].edges.getByBoundingBox(
        xMin=-tol, yMin=-tol, zMin=-tol, xMax=tol, yMax=Hr*voxY+tol, zMax=tol) , 
        size=voxY, constraint=FIXED)
    del mdb.models[modelo].sketches['__profile__']
    mdb.models[modelo].parts['block'].seedEdgeBySize(
        edges=mdb.models[modelo].parts['block'].edges.getByBoundingBox(
        xMin=-tol, yMin=Hr*voxY-tol, zMin=-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=tol) , 
        size=voxX, constraint=FIXED)
    mdb.models[modelo].parts['block'].seedEdgeBySize(
        edges=mdb.models[modelo].parts['block'].edges.getByBoundingBox(
        xMin=Wr*voxX-tol, yMin=-tol, zMin=-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=tol) , 
        size=voxY, constraint=FIXED)
    mdb.models[modelo].parts['block'].seedEdgeBySize(
        edges=mdb.models[modelo].parts['block'].edges.getByBoundingBox(
        xMin=-tol, yMin=-tol, zMin=-tol, xMax=tol, yMax=tol, zMax=Lr*voxZ+tol) , 
        size=voxZ, constraint=FIXED)
    #mdb.models[modelo].parts['block'].setMeshControls(elemShape=QUAD, regions=
    #    mdb.models[modelo].parts['block'].cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
    #elemType1 = mesh.ElemType(elemCode=CPS4R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=DEFAULT, elemDeletion=ON)  ### BM FLAG
    elemType1 = mesh.ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    pickedRegions =(mdb.models[modelo].parts['block'].cells.getSequenceFromMask(mask=('[#1 ]', ), ), )
    mdb.models[modelo].parts['block'].setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
    mdb.models[modelo].parts['block'].generateMesh()
    mdb.models[modelo].parts['block'].Set(name='set-check',elements = mdb.models[modelo].parts['block'].elements[0:int(math.floor(Hr*Wr/3.+Lr*5.))])
    mdb.models[modelo].rootAssembly.Instance(name='Lattice', part=mdb.models[modelo].parts['block'], dependent=ON)
    
    print 'DO YOU MOVE IN A KIND OF SLEEP?'
    keepnodes = []
    #label_rho = np.zeros(sum(propMatC))
    for mm in range(0,int(tot_props)):
        if propMatC[mm]>0.0:
            indlabels=materialMat[mm,0:int(propMatC[mm])].tolist()
            denss = prop_list[mm]
            if denss == 0.01:
                denss = 0.0
            #label_rho[(np.array(indlabels)-1)]= denss
            elems=mdb.models[modelo].parts['block'].elements.sequenceFromLabels(indlabels)
            #if mm == 2.0:
            mdb.models[modelo].parts['block'].Set(elements=elems,name='elems-'+str(mm))  #this group will be empty by the second iteration!
            
            eps_ult =  (fit_eps_ult[0]*np.exp(-denss*fit_eps_ult[1]) + fit_eps_ult[2]  )*punish
            # A, B, alpha, k, beta
            #A_a = fit_A[0]*denss**fit_A[1] + fit_A[2]
            B_a =    1.0
            if denss < 0.4 :
                rho_design = denss**np.transpose(np.array(([range(0,order_A+1)])))
                A_a =    np.dot(np.transpose(rho_design),theta_A)[0,0]
                A_a40 = fit_A40_1[0]*denss**fit_A40_1[1] + fit_A40_1[2]
                alpha_a40 = fit_alpha40_1[0]*denss**fit_alpha40_1[1] + fit_alpha40_1[2]
            else:
                A_a = fit_A2[0]*denss**fit_A2[1] + fit_A2[2]
                A_a40 = fit_A40_2[0]*denss**fit_A40_2[1] + fit_A40_2[2]
                alpha_a40 = fit_alpha40_2[0]*denss**fit_alpha40_2[1] + fit_alpha40_2[2]
            
            if denss < 0.2 :
                alpha_a = fit_alpha1[0]*denss**fit_alpha1[1] + fit_alpha1[2]
                A_a20 = fit_A20_1[0]*denss**fit_A20_1[1] + fit_A20_1[2]
                alpha_a20 = fit_alpha20_1[0]*denss**fit_alpha20_1[1] + fit_alpha20_1[2]
            else:
                alpha_a = fit_alpha2[0]*denss**fit_alpha2[1] + fit_alpha2[2]
                A_a20 = fit_A20_2[0]*denss**fit_A20_2[1] + fit_A20_2[2]
                alpha_a20 = fit_alpha20_2[0]*denss**fit_alpha20_2[1] + fit_alpha20_2[2]
					
            A_a =  A_a40
            alpha_a =  alpha_a40
            K_a = fit_ratio[0]*np.exp(-denss*fit_ratio[1]) + fit_ratio[2]
            beta_a = alpha_a/K_a
            #initial Ym
            Ym_liu = A_a*alpha_a/2.
            if denss<0.2:
                young = Ym_liu
            elif denss>=0.2:
                young = Ym_liu
            ##else:
            #
            # define the plastic strain vector
            eps_off = 0.02
            eps = np.transpose(np.array(([np.linspace(0,eps_ult+eps_off,pointsnum)])))
            sig = A_a*( np.exp(alpha_a*eps) - 1. )/( B_a + np.exp(beta_a*eps) )
            # plastic stress vector
            eps_pl = eps - sig/young -eps_off
            if len( np.where(eps_pl<0)[0] ) != 0 :
                ind_incorr = np.where(eps_pl<0)[0][-1] 
                eps_pl = eps_pl[ind_incorr::] - eps_pl[ind_incorr]
                eps = eps[ind_incorr::]
                sig_pl = sig[ind_incorr::]
                sig = sig[ind_incorr::] 
                sig_y_a = sig_pl[0][0] 
            # Add extra points if there's no plastic section
            eps_kill = 0.05
            if len( eps_pl) == 1 :
                eps_pl = np.concatenate((eps_pl, np.array(([ [  punish*eps_kill+eps_pl[-1][0]  ]]))))
                young = sig_y_a/eps_ult
                sig_pl = np.concatenate((sig_pl, np.array(([ [  sig_y_a + 0.99*young*eps_pl[-1][0]  ]]))))
            
            if eps_pl.size == 0:
                sig_pl = np.array(([sig[-1,:]]))
                eps_pl = np.array(([0.0]))
            else:
                eps_pl[0] = 0.0
            killlllls = np.zeros(eps_pl.size)
            for zzz in range(1,len(eps_pl)): 
                if eps_pl[zzz] < eps_pl[zzz-1]: 
                    eps_pl[zzz] = eps_pl[zzz-1]
                    killlllls[zzz] = 1
            
            corr_inds = (np.where(killlllls==0))[0]
            sig_pl = sig_pl[killlllls==0]
            eps_pl = eps_pl[killlllls==0]
            k_l =  voxY
            eps_pl_D_ini = (eps_pl[-1] - eps_kill*punish)[0]
            if eps_pl_D_ini <= 0  :
                eps_pl_D_ini = 0.001
                eps_kill = (eps_pl[-1])[0]
            u_pl = np.concatenate((np.array(([ [  0.0  ]])), np.array(([ [  punish*eps_kill*k_l  ]]))))
            D_vec = D_fin*u_pl/u_pl.max()
            ##yeah and that's it
            #oopsi
            if mm > 0:
                elset = mdb.models[modelo].parts['block'].sets['elems-'+str(mm)]
                elset_nodes = set()
                for element in elset.elements:
                    elset_nodes.update(element.connectivity)
                keepnodes = keepnodes + list(elset_nodes)
                mdb.models[modelo].Material(name='Mat-'+str(mm))
                mdb.models[modelo].materials['Mat-'+str(mm)].Density(table=((denss, ), ))
                poissons = denss*(0.4) + (1-denss)*0.49
                mdb.models[modelo].HomogeneousSolidSection(material='Mat-'+str(mm), name='Sec-'+str(mm),thickness=None)#*0.03/0.027)
                if crackBC == 0 and denss < 0.05 : # Soft no crack is hyperelastic material!
                    mdb.models[modelo].materials['Mat-'+str(mm)].Hyperelastic( materialType=ISOTROPIC, testData=OFF, type=OGDEN, volumetricResponse=VOLUMETRIC_DATA, table=((mu_soft_exp, alpha_soft_exp, D1), ))
                    #arr = np.transpose(((sig_pl_soft_HE[:]),(eps_pl_soft_HE[:])))
                    #mdb.models[modelo].materials['Mat-'+str(mm)].Plastic(table=tuple(map(tuple, arr)))
                if crackBC == 0 and denss >= 0.05 : # Only plastic for the rest materials
                    mdb.models[modelo].materials['Mat-'+str(mm)].Elastic(table=((young, poissons), ))
                    arr = np.single(np.transpose((   (np.squeeze(sig_pl[:])) , (np.squeeze(eps_pl[:]))   )))
                    mdb.models[modelo].materials['Mat-'+str(mm)].Plastic(table=tuple(map(tuple, arr)))
                if crackBC == 1 : # Soft crack is elastic but with the final true value material!
                    mdb.models[modelo].materials['Mat-'+str(mm)].Elastic(table=((young, poissons), ))
                    arr = np.single(np.transpose((   (np.squeeze(sig_pl[:])) , (np.squeeze(eps_pl[:]))   )))
                    mdb.models[modelo].materials['Mat-'+str(mm)].Plastic(table=tuple(map(tuple, arr)))
                    mdb.models[modelo].materials['Mat-'+str(mm)].DuctileDamageInitiation(    table=((eps_pl_D_ini, -0.3333, 0.0), ))
                    arr = np.single(np.transpose((   (np.squeeze(D_vec[:])) , (np.squeeze(u_pl[:]))   )))
                    mdb.models[modelo].materials['Mat-'+str(mm)].ductileDamageInitiation.DamageEvolution(    type=DISPLACEMENT, softening=TABULAR, table=tuple(map(tuple, arr))  )
    
    #### NOW WE CREATE PROPERTIES: AN ARRAY THAT HAS THE TIMES EACH OF THE N*M+1 PROPERTIES APPEAR
    print 'TIME HAS SLIPPED AWAY.'
    nodenum = len(mdb.models[modelo].parts['block'].nodes)
    nodlabs = np.arange(1,nodenum+1)
    knlabs = np.sort(np.asarray(keepnodes))
    killlabs = np.delete(nodlabs,knlabs)#-1)  
    nodess=mdb.models[modelo].parts['block'].nodes.sequenceFromLabels(killlabs)
    mdb.models[modelo].parts['block'].Set(nodes=nodess,name='nodeskill')
    mdb.models[modelo].parts['block'].PartFromMesh(copySets=True, name='block-mesh')
    mdb.models[modelo].parts['block-mesh'].deleteElement(elements=mdb.models[modelo].parts['block-mesh'].sets['elems-0'])
    mdb.models[modelo].parts['block-mesh'].deleteNode(nodes=mdb.models[modelo].parts['block-mesh'].sets['nodeskill'])
    del mdb.models[modelo].parts['block-mesh'].sets['elems-0']
    del mdb.models[modelo].parts['block-mesh'].sets['nodeskill']
    del keepnodes, nodlabs, nodenum, knlabs, killlabs, elset_nodes, indlabels
    #oopsi
    mdb.models[modelo].rootAssembly.Instance(name='Lattice', part=mdb.models[modelo].parts['block-mesh'], dependent=ON)
    mdb.models[modelo].rootAssembly.Set(name='ele',elements=[mdb.models[modelo].rootAssembly.instances['Lattice'].elements.getByBoundingBox(xMin=0-tol, yMin=0-tol, zMin=0.0-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=voxZ*Lr+tol)],)
    mdb.models[modelo].rootAssembly.Set(name='eleAll',elements=[mdb.models[modelo].rootAssembly.instances['Lattice'].elements.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=voxZ*Lr+tol)],)
    ele_len = len(mdb.models[modelo].rootAssembly.sets['ele'].elements)
    eleAll_len = len(mdb.models[modelo].rootAssembly.sets['eleAll'].elements)
    #lab2ele = np.zeros(ele_len)
    ele2lab = np.zeros(ele_len)
    for ee in range(0,ele_len):
        #lab2ele[mdb.models[modelo].rootAssembly.sets['ele'].elements[ee].label] = ee
        ele2lab[ee] = mdb.models[modelo].rootAssembly.sets['ele'].elements[ee].label
    
    a=mdb.models[modelo].parts['block-mesh']
    a.regenerate()
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(    meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(    referenceRepresentation=OFF)
    session.viewports['Viewport: 1'].enableMultipleColors()
    session.viewports['Viewport: 1'].setColor(initialColor='#BDBDBD')
    for mm in range(1,int(tot_props)):
      if propMatC[mm]>0.0:
        indlabels=materialMat[mm,0:int(propMatC[mm])].tolist()
        elems=mdb.models[modelo].parts['block-mesh'].elements.sequenceFromLabels(indlabels)
        mdb.models[modelo].parts['block-mesh'].Set(elements=elems,name='elems-'+str(mm))  #this group will be empty by the second iteration!   
        mdb.models[modelo].parts['block-mesh'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(elements=elems), 
            sectionName='Sec-'+str(mm),thicknessAssignment=FROM_SECTION)
        denss = prop_list[mm]
        ##hexstr = (str(hex(      255-int(round(denss*255))    )[-2:]).upper()+str(hex(   00    )[-2:]).upper()+str(hex(int(round(    denss*255    )))[-2:]).upper()).replace('X', '0')
        hexstr = (str(hex(      249-int(round(denss*(249-0)))    )[-2:]).upper()+str(hex(   249-int(round(denss*(249-102)))    )[-2:]).upper()+str(hex(int(round(    255-int(round(denss*(255-162)))     )))[-2:]).upper()).replace('X', '0')
        cmap = session.viewports['Viewport: 1'].colorMappings['Section']
        cmap.updateOverrides(overrides={'Sec-'+str(mm):(True, '#FF0000', 'Default', '#'+hexstr)})
    
    print 'YOUR LIFE IS STOLEN.'
    #session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(meshVisibleEdges=FREE)
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
    session.viewports['Viewport: 1'].setColor(colorMapping=cmap)
    session.viewports['Viewport: 1'].disableMultipleColors()
    session.graphicsOptions.setValues(backgroundStyle=SOLID,     backgroundColor='#000000')
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(    meshVisibleEdges=FREE)
    #session.printToFile(fileName=picfolder + '\\D-' + modelo , format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
    ##################################################################################
    #FINISHED CREATING THE DEVICE, TIME TO CREATE STEP AND APPLY BC
    ##################################################################################
    mdb.models[modelo].rootAssembly.Set(name='Top',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=voxZ*Lr-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=voxZ*Lr+tol),],)
    mdb.models[modelo].rootAssembly.Set(name='Bot',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=0.0+tol),],)
    mdb.models[modelo].rootAssembly.Set(name='Left',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=0.0+tol, yMax=Hr*voxY+tol, zMax=voxZ*Lr+tol),],)
    mdb.models[modelo].rootAssembly.Set(name='Right',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=Wr*voxX-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=Wr*voxX+tol, yMax=Hr*voxY+tol, zMax=voxZ*Lr+tol),],)
    mdb.models[modelo].rootAssembly.Set(name='Ymirr',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=Wr*voxX+tol, yMax=0.0+tol, zMax=voxZ*Lr+tol),],)
    
    #### CREATE THE STEP AND ADD THE BCs
    print 'YOU TARRIED WITH TRIFFLES.'
    ##################################################################################
    #CREATE STEP AND APPLY BC
    ##################################################################################
    mdb.models[modelo].StaticStep(name='Step-1', nlgeom=ON, previous='Initial',maxNumInc=666, initialInc=0.1,minInc=1E-8,maxInc=1.0,)
    mdb.models[modelo].steps['Step-1'].setValues( matrixStorage=UNSYMMETRIC )
    mdb.models[modelo].steps['Step-1'].control.setValues(allowPropagation=OFF, resetDefaultValues=OFF, timeIncrementation=(4.0, 10.0, 9.0, 16.0, 10.0, 4.0, 12.0, 10.0, 6.0, 3.0, 50.0))
    mdb.models[modelo].steps['Step-1'].control.setValues(discontinuous=ON, timeIncrementation=(8.0, 10.0, 9.0, 16.0, 10.0, 4.0, 12.0, 12.0, 6.0, 3.0, 50.0))
    alpha = strainperc/100.*Wr*voxX
    mdb.models[modelo].rootAssembly.Set(name='PerBound',
        nodes=[
        mdb.models[modelo].rootAssembly.instances['Lattice'].nodes.getByBoundingBox(xMin=0.0-tol, yMin=0.0-tol, zMin=0.0-tol, xMax=0.0+tol, yMax=0.0+tol, zMax=0.0+tol),
        ],)
    (CoorFixNode,NameRef1, NameRef2, NameRef3)=PeriodicBound3D_reduced(mdb,modelo,'PerBound',[(Wr*voxX,0.0,0.0),(0.0,Hr*voxY,0.0),(0.0,0.0,voxZ*Lr)],[0,1,1])
    mdb.models[modelo].rootAssembly.translate(instanceList=('RefPoint-1', ), vector=(Wr*voxX, Hr*voxY*0.5, Lr*voxZ*0.5))
    regionDef=mdb.models[modelo].rootAssembly.sets['RefPoint-1']
    mdb.models[modelo].historyOutputRequests['H-Output-1'].setValues(variables=('U1','U2','U3','RF1','RF2','RF3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE, frequency=1,)
    mdb.models[modelo].fieldOutputRequests['F-Output-1'].setValues( variables=('S','TRIAX','PE','PEEQ','PEMAG','LE','U','RF','DAMAGEC','DAMAGET','DAMAGESHR','SDEG','DMICRT','SHRRATIO','STATUS'))
    region1=mdb.models[modelo].rootAssembly.sets[NameRef2]
    region2=mdb.models[modelo].rootAssembly.sets['Right']
    mdb.models[modelo].Coupling(name='Constraint-RN', controlPoint=region1, surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)  
    mdb.models[modelo].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name= 'BC-REF-0', region=mdb.models[modelo].rootAssembly.sets['Left'],  u1=0., u2=0., u3=0., ur1=0.0, ur2=0.0, ur3=0.0)
    mdb.models[modelo].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name= 'BC-REF-1', region=mdb.models[modelo].rootAssembly.sets['RefPoint-1'],  u1=alpha, u2=0.0, u3=0., ur1=0.0,ur2=0.0,ur3=0.0)
    
    ##################################################################################
    #JOB AND RUN
    ##################################################################################
    
    #### CREATE THE JOB AND ADD IT TO THE BATCH FILE, THEN DO THE ABBA FILES
    print 'VICTIM OF YOUR FOLLY.'
    mdb.Job(name=modelo, model=modelo, description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=CPUSS, 
        numDomains=6, numGPUs=0)
    
    mdb.jobs[modelo].writeInput(consistencyChecking=OFF) 

session.Spectrum(name="Spectrum-1", colors =('#7C10EA', '#7653FB', '#668FFB',  '#51BFEC', '#53E5C6', '#75FA9B', '#9DFF5C', '#CBF41E', '#EBD701',  '#FEAA00', '#FB7501', '#EA3300', ))
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(  outsideLimitsMode=SPECIFY, outsideLimitsAboveColor='#930E00',   outsideLimitsBelowColor='#7D02DE')
session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(spectrum='Spectrum-1')
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(    visibleEdges=FREE)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])


print ' '
print '    THOU SHALT NOT MAKE A MACHINE IN THE LIKENESS OF A MAN''S MIND.'
print ' '


