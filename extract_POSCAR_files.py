#--------
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#------------------------------

for k in range(2):
   if (k == 0):
      database_file = 'SAMBA_Monolayers_database.json'
      POSCAR_files = 'POSCAR_Monolayers'
   if (k == 1):
      database_file = 'SAMBA_Bilayers_database.json'
      POSCAR_files = 'POSCAR_Bilayers'

   #------------------------------------------------------------
   with open(database_file, "r") as file: data = json.load(file)        # Carregando o banco de dados
   df = pd.json_normalize(data if isinstance(data, list) else [data])   # Convertendo os dados para um DataFrame do pandas
   tags = df.columns                                                    # Identificando as diferentes classes de tags
   #------------------------------------------------------------
   if (k == 1):
      print(f'Número total de entradas no banco de dados: {len(data)}')
      print(f'Número total de diferentes tags (chaves): {len(tags)}')
      print(" ")
      print("Classes de Tags no banco de dados:")
      print(tags)
      #----------

   #--------------------------------------
   if os.path.isdir(POSCAR_files): 0 == 0
   else: os.mkdir(POSCAR_files)
   #-----------------------------
   for material in data:
       #-------------------------
       label = material.get("id")   
       vetor_A1 = material.get("a1")
       vetor_A2 = material.get("a2")
       vetor_A3 = material.get("a3") 
       type_ions = material.get("type_ions_layers") 
       type_nions = material.get("number_type_ions_layers")
       direct_coord = material.get("direct_coord_ions")
       #------------------------------------------------------------
       poscar = open(POSCAR_files + '/' + str(label) + '.vasp', "w")
       poscar.write(f'{str(label)} \n')
       poscar.write(f'1.0 \n')
       poscar.write(f'{vetor_A1[0]} {vetor_A1[1]} {vetor_A1[2]} \n')
       poscar.write(f'{vetor_A2[0]} {vetor_A2[1]} {vetor_A2[2]} \n')
       poscar.write(f'{vetor_A3[0]} {vetor_A3[1]} {vetor_A3[2]} \n')
       #------------------------------
       for i in range(len(type_ions)):
           for j in range(len(type_ions[i])):
               poscar.write(f'{type_ions[i][j]} ')
       poscar.write(f' \n')
       #-------------------------------
       for i in range(len(type_nions)):
           for j in range(len(type_nions[i])):
               poscar.write(f'{type_nions[i][j]} ')
       poscar.write(f' \n')
       #-------------------------
       poscar.write(f'Direct \n')
       #---------------------------------
       for i in range(len(direct_coord)): poscar.write(f'{direct_coord[i][0]} {direct_coord[i][1]} {direct_coord[i][2]} \n')
       #---------------------------------
       poscar.close()
       #-------------

print(" ")
print("==================================")
print("Arquivos POSCAR extraídos: =======")
print("==================================")
