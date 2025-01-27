import rmsd as rmsdlib
from pathlib import Path
import numpy as np
import tempfile
import os


def read_chunks(filename, chunk_size=17, reference_files=None):
   """
   Read and process XYZ file in chunks and collect similar structures
   """
   # Initialize counters and output files
   structure_counts = {
       'a': 0,
       'b': 0,
       'c': 0,
       'd': 0
   }
   
   # Open output files
   output_files = {
       'a': open('a.xyz', 'w'),
       'b': open('b.xyz', 'w'), 
       'c': open('c.xyz', 'w'),
       'd': open('d.xyz', 'w')
   }
   
   if reference_files is None:
       reference_files = {
           'a': 'ci_a.xyz',
           'b': 'ci_b.xyz',
           'c': 'ci_c.xyz',
           'd': 'ci_d.xyz'
       }

   try:
       with open(filename, 'r') as file:
           chunk_number = 0
           while True:
               chunk = []
               for _ in range(chunk_size):
                   line = file.readline()
                   if not line:
                       break
                   chunk.append(line)
               
               if not chunk:
                   break
               
               chunk_number += 1
               with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.xyz') as temp_file:
                   temp_file.writelines(chunk)
                   temp_path = temp_file.name

               try:
                   filenames = {**reference_files, 'e': temp_path}
                   atoms = {}
                   coord = {}
                   
                   for letter in ['a', 'b', 'c', 'd', 'e']:
                       atoms[letter], coord[letter] = rmsdlib.get_coordinates_xyz(filenames[letter])
                       coord[letter] -= rmsdlib.centroid(coord[letter])

                   rmsd_method = rmsdlib.kabsch_rmsd
                   reorder_method = rmsdlib.reorder_hungarian
                   rmsd = {}
                   
                   for letter in ['a', 'b', 'c', 'd']:
                       rmsd[letter], _, _, _ = rmsdlib.check_reflections(
                           atoms[letter],
                           atoms['e'],
                           coord[letter],
                           coord['e'],
                           reorder_method=reorder_method,
                           rmsd_method=rmsd_method,
                       )
                   
                   min_key = min(rmsd, key=rmsd.get)
                   structure_counts[min_key] += 1
                   
                   # Write to appropriate output file
                   output_files[min_key].writelines(chunk)
                   print(f"Chunk {chunk_number}: Best match = {min_key} (RMSD = {rmsd[min_key]:.4f})")
               
               finally:
                   os.unlink(temp_path)
                   
       # Print final counts
       print("\nFinal Structure Counts:")
       print("-" * 20)
       total = sum(structure_counts.values())
       for structure, count in structure_counts.items():
           percentage = (count / total) * 100 if total > 0 else 0
           print(f"Structure {structure}: {count} ({percentage:.1f}%)")
       print(f"Total geometries analyzed: {total}")
           
   except FileNotFoundError:
       print(f"Error: File '{filename}' not found")
   except Exception as e:
       print(f"Error reading file: {e}")
   finally:
       # Close all output files
       for f in output_files.values():
           f.close()
       
   return structure_counts

# Example usage
filename = "pre_hop_all.xyz"
counts = read_chunks(filename)