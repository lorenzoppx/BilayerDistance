from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import JmolNN
#---------------------------------------------------------------
import numpy as np
import filecmp
import shutil
import sys  
import os

#------------------------
fator_supercell = [5,5,5]
raio_corte = 10.0
#----------------

#--------------------------------------------------------------------------------------------------------------------
# Rotulo temporário para os ions inequivalentes da célula unitária --------------------------------------------------
alphabet  = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','X','Z','Y','W']
alphabet += ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','x','z','y','w'] 
#-------------------------------------------------------------------------------------------------------------------- 

for ii in range(2):
    #----------------------------------------
    if (ii == 0): diret = 'POSCAR_Monolayers'
    if (ii == 1): diret = 'POSCAR_Bilayers'
    #------------------------------------------------
    # Listando todos os arquivos com extensão '.vasp'
    poscar_files = []
    for root, _, files in os.walk(diret):
        for file in files:
            if file.endswith('.vasp'): poscar_files.append(os.path.join(file))
    #-------------------------------------------------------------------------
    for iii in range(len(poscar_files)):
        #----------------------------------------
        Lattice = diret + '/' + poscar_files[iii]
        print(f'Analisando: {Lattice}')
        #------------------------------

        #------------------------------------------------
        ion_label  = []; ion_label_string  = ''
        ion_label_temp  = []; ion_label_temp_string  = ''
        nlabel = []; nions = 0
        #------------------------------------------------



        #===============================================================================
        # Reescrita do arquivo POSCAR, de forma a separar os diferentes ions da rede ===
        #===============================================================================
        poscar = open(Lattice, "r")
        #--------------------------
        for i in range(5):
            VTemp = poscar.readline()
        VTemp1 = poscar.readline().split()
        VTemp2 = poscar.readline().split()
        #---------------------------------
        STemp1 = ''
        STemp2 = ''
        number2 = 0
        for i in range(len(VTemp2)):
            number = 0 
            for j in range(int(VTemp2[i])):
                number += 1
                number2 += 1
                STemp1 += str(VTemp1[i]) + '_' + str(number2) + ' '
                STemp2 += '1 '
        #---------------------
        poscar.close() 
        #-------------

        with open(Lattice, 'r') as arquivo: lines = arquivo.readlines()
        lines[6 -1] = STemp1 + '\n'
        lines[7 -1] = STemp2 + '\n'
        with open('POSCAR_temp1.vasp', 'w') as arquivo: arquivo.writelines(lines)



        #===============================================================================
        # Reescrita do arquivo POSCAR, de forma a separar os diferentes ions da rede ===
        #===============================================================================
        poscar = open('POSCAR_temp1.vasp', "r")
        poscar_new = open('POSCAR_temp2.vasp', "w")

        for i in range(5):
            VTemp = poscar.readline()
            poscar_new.write(f'{VTemp}')

        #=================================================
        VTemp = poscar.readline().split()
        for i in range(len(VTemp)):  
            ion_label.append(str(VTemp[i]))                    # Armazenando o label de cada ion da rede (individualmente) na posição de um vetor.        
            ion_label_string += str(VTemp[i]) + ' '            # Criando uma String única contendo os labels de todos os ions da rede.         
        ion_label_string = ion_label_string[:-1]
        #-------------------------------------------------
        VTemp = poscar.readline().split()
        for i in range(len(VTemp)):                         
            nlabel.append(int(VTemp[i]))                       # Armazenando o número de cada tipo de ion da rede (individualmente) na posição de um vetor.
            nions += int(VTemp[i])                             # Obtendo o nº total de ions da rede.
        #-------------------------------------------------
        cont = -1;
        for i in range(len(ion_label)):                        # Loop sobre os diferentes tipos de ions da rede.
            ion_label_temp.append('')                          # Para cada diferente tipo de ion da rede, cria-se uma String.
            temp_label = ''
            for j in range(nlabel[i]):                         # Loop sobre o número de cada tipo de ion da rede.
                cont += 1
                temp_label += alphabet[cont] + ' '                                           
            ion_label_temp_string += temp_label        
            ion_label_temp[i] += temp_label[:-1]               # Criando uma String (para cada ion da rede individualmente) que armazena os rotulos temporários das correspondentes subredes.
        ion_label_temp_string = ion_label_temp_string[:-1]     # Criando uma String única que armazena os rotulos de todas as subredes.
        #-------------------------------------------------
        vector_ions_labels    = ion_label                 # Armazenando o label de cada ion da rede (individualmente) na posição de um vetor.
        string_ions_labels    = ion_label_string          # Armazenando a String contendo os labels de todos os ions da rede.
        vector_subrede_labels = ion_label_temp            # Armazenando os labels de subrede referentes a cada ion da rede (individualmente) na posição de um vetor.
        string_subrede_labels = ion_label_temp_string     # Armazenando a String contendo os labels de subrede referentes a todos os ions da rede.
        #-------------------------------------------------
        for i in range(nions):
            poscar_new.write(f'{alphabet[i]} ')
        poscar_new.write("\n")
        #-------------------------------------------------
        for i in range(nions):
            poscar_new.write("1 ")
        poscar_new.write("\n")
        #-------------------------------------------------
        VTemp = poscar.readline()
        poscar_new.write(f'{VTemp}')
        #-------------------------------------------------
        for i in range(nions+1):
            VTemp = poscar.readline()
            poscar_new.write(f'{VTemp}')
        #-------------------------------------------------
        poscar.close()
        poscar_new.close()
        #=================================================
        os.remove('POSCAR_temp1.vasp')
        #-----------------------------



        #==========================================================================
        # Criando a Supercélula das redes (em coordenadas diretas) ================
        #==========================================================================
        structure = Poscar.from_file('POSCAR_temp2.vasp').structure
        # Cria uma supercélula multiplicando os vetores da rede
        supercell = structure.copy()
        supercell.make_supercell(fator_supercell)
        Poscar(supercell).write_file('POSCAR_Supercell_direto.vasp')
        #=====================================================================================
        # Evitando erro na escrita do rotulo das subredes (bug no pymatgen) ==================
        #=====================================================================================
        with open('POSCAR_Supercell_direto.vasp', 'r') as arquivo: lines = arquivo.readlines()
        lines[6 -1] = string_subrede_labels + '\n'
        with open('POSCAR_Supercell_direto.vasp', 'w') as arquivo: arquivo.writelines(lines)
        #=====================================================================================
        os.remove('POSCAR_temp2.vasp')
        #-----------------------------



        #==========================================================================
        # Convertendo as coordenadas da Supercélula para a forma cartesiana =======
        #==========================================================================

        #------------------------------------------------------
        # Obtenção dos novos vetores da rede ------------------
        #------------------------------------------------------
        poscar = open('POSCAR_Supercell_direto.vasp', "r")
        VTemp = poscar.readline()
        VTemp = poscar.readline(); param = float(VTemp)
        VTemp = poscar.readline().split(); A = [float(VTemp[0])*param, float(VTemp[1])*param, float(VTemp[2])*param]
        VTemp = poscar.readline().split(); B = [float(VTemp[0])*param, float(VTemp[1])*param, float(VTemp[2])*param]
        VTemp = poscar.readline().split(); C = [float(VTemp[0])*param, float(VTemp[1])*param, float(VTemp[2])*param]
   
        #----------------------------------------------------------------------
        # Obtenção das coordenadas da nova origem do sistema de coordenadas ---
        #----------------------------------------------------------------------
        center_x = 0.5*A[0] + 0.5*B[0] + 0.5*C[0]
        center_y = 0.5*A[1] + 0.5*B[1] + 0.5*C[1]
        center_z = 0.5*A[2] + 0.5*B[2] + 0.5*C[2]

        #-----------------------------------------------------------------
        # Armazenando os rótulos das subredes (em um vetor de strings) ---
        #-----------------------------------------------------------------
        VTemp = poscar.readline().split()
        vector_ion = []
        for i in range(len(VTemp)):
            vector_ion.append(str(VTemp[i]))
        #-----------------------------------------------------------------------------------------
        # Armazenando o nº de ions de cada subrede (em um vetor), e o nº total de ions da rede ---
        #-----------------------------------------------------------------------------------------
        VTemp = poscar.readline().split()
        vector_nion = []; passo = 0
        for i in range(len(VTemp)):
            vector_nion.append(str(VTemp[i]))
            passo += int(VTemp[i])
        #-------------
        poscar.close()
        #-----------------------------------------------------------------------------------
        # Criando um vetor que associa a cada ion da rede o correspondente rotulo de subrede
        #-----------------------------------------------------------------------------------
        vector_rot_subredes = ['0']*passo
        number = -1
        for i in range(len(vector_nion)):
            for j in range(int(vector_nion[i])):
                number += 1 
                vector_rot_subredes[number] = str(vector_ion[i])
                #-----------------------------------------------

        #---------------------------------------------------------
        # Escrita do arquivo POSCAR em coordenadas cartesianas ---
        #---------------------------------------------------------
        poscar = open('POSCAR_Supercell_direto.vasp', "r")
        poscar_new = open('POSCAR_Supercell_cartesiano.vasp', "w")
        #---------------------------------------------------------
        poscar_new.write(f'{ion_label_string} \n')
        VTemp = poscar.readline()
        #------------------------
        for i in range(6):
            VTemp = poscar.readline()
            poscar_new.write(f'{VTemp}')
        VTemp = poscar.readline()
        poscar_new.write(f'cartesian \n')
        #--------------------------------
        for j in range(passo):
            VTemp = poscar.readline().split()
            k1 = float(VTemp[0]); k2 = float(VTemp[1]); k3 = float(VTemp[2])
            coord_x = ((k1*A[0]) + (k2*B[0]) + (k3*C[0]))*param
            coord_y = ((k1*A[1]) + (k2*B[1]) + (k3*C[1]))*param
            coord_z = ((k1*A[2]) + (k2*B[2]) + (k3*C[2]))*param
            poscar_new.write(f'{coord_x:>28,.21f} {coord_y:>28,.21f} {coord_z:>28,.21f}  {vector_rot_subredes[j]} \n')
        #-----------------
        poscar.close()
        poscar_new.close()
        #-----------------
        os.remove('POSCAR_Supercell_direto.vasp')
        #----------------------------------------



        #-------------------------------------------------------------
        poscar_file = poscar_files[iii].replace('.vasp', '_coord.txt')
        output = open(diret + '/' + poscar_file, "w")
        #--------------------------------------------

        for k in range(nions):
            #---------------------------
            tipo_ion = ion_label_temp[k]

            #====================================================
            # Obtendo o nº de coordenação de cada ion da rede ===
            #====================================================
            structure = Structure.from_file(Lattice)                           # Carregar a estrutura a partir do arquivo POSCAR
            sg = StructureGraph.with_local_env_strategy(structure, JmolNN())   # Criando o grafo da estrutura com base na estratégia JmolNN
            #---------------------------------------------------------------                                 
            coordination = len(sg.get_connected_sites(k))                      # Obtendo a coordenação do sítio
            output.write(f'{ion_label[k]} = [ {coordination}, ')



            #============================================================================
            # Encontrando o ion mais próximo da nova origem do sistema de coordenadas ===
            #============================================================================

            #-------------------------------------------------------------------------
            # Obtenção das coordenadas do átomo mais próximo da nova origem ----------
            #-------------------------------------------------------------------------
            poscar = open('POSCAR_Supercell_cartesiano.vasp', "r")
            for i in range(8): VTemp = poscar.readline()
            temp_d = 1000.0
            for i in range(passo):
                VTemp = poscar.readline().split()
                #------------------------
                coord_x = float(VTemp[0])
                coord_y = float(VTemp[1])
                coord_z = float(VTemp[2])
                ion     = str(VTemp[3]) 
                #------------------------
                dist = ((coord_x -center_x)**2 + (coord_y -center_y)**2 + (coord_z -center_z)**2)**(0.5)
                if ((dist <= temp_d) and (ion == tipo_ion)):
                   temp_d = dist 
                   new_center_x = coord_x; new_center_y = coord_y; new_center_z = coord_z
            #-------------
            poscar.close()
            #-------------

            #--------------------------------------------------------------------------
            # Deslocando as coordenadas cartesianas e Armazenando as coordenações -----
            #--------------------------------------------------------------------------
            poscar = open('POSCAR_Supercell_cartesiano.vasp', "r")
            poscar_new = open('Coordenacao_' + str(tipo_ion) + '.vasp', "w")
            #---------------------------------------------------------------
            for i in range(8):  VTemp = poscar.readline()
            for i in range(passo):
                VTemp = poscar.readline().split()
                #---------------------------------------
                coord_x = float(VTemp[0]) - new_center_x
                coord_y = float(VTemp[1]) - new_center_y
                coord_z = float(VTemp[2]) - new_center_z
                label_ion = str(VTemp[3])
                #-----------------------------------------------------------------------
                dist = ((coord_x)**2 + (coord_y)**2 + (coord_z)**2)**(0.5)
                if (dist <= raio_corte): poscar_new.write(f'{dist:>10,.6f}_{label_ion[0]} \n')
            #-----------------
            poscar.close()
            poscar_new.close()
            #-----------------

            file = np.loadtxt('Coordenacao_' + str(tipo_ion) + '.vasp', dtype='str') 
            file.shape

            x = file
            coord, n_coord = np.unique(x, return_counts=True)

            ncoord = 0 
 
            for j in range(len(coord)):
                #----------------------------------------
                VTemp = coord[j].replace('_',' ').split()
                #----------------------------------------
                for m in range(len(ion_label)):
                    if (tipo_ion == ion_label_temp[m]): temp1 = ion_label[m]
                    if (VTemp[1] == ion_label_temp[m]): temp2 = ion_label[m]  

                if (VTemp[0] != '0.000000'):
                   ncoord += int(n_coord[j])
                   output.write(f'[{n_coord[j]}, {temp2}, {VTemp[0]}], ')
                   # if (ncoord <  2*coordination): output.write(f', ')                  
            output.write(f' ] \n')

            #--------------------------------------------------
            os.remove('Coordenacao_' + str(tipo_ion) + '.vasp')
        output.close() 
        os.remove('POSCAR_Supercell_cartesiano.vasp')
        #--------------------------------------------
