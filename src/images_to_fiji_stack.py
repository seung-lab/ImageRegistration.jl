# images_to_fiji_stack.py
# Load image directory into FIJI, apply translations, and merge into stack

import os
import csv

resize_factor = 0.25

bucket_dir_path = '/usr/people/tmacrina/seungmount'
datasets_dir_path = "research/Julimaps/datasets"
cur_dataset = "piriform"

premontaged_dir_path = "1_premontaged"
montaged_dir_path = "2_montaged"
prealigned_dir_path = "3_prealigned"
aligned_dir_path = "4_aligned"

premontaged_offsets_filename = 'premontaged_offsets.txt'
montaged_offsets_filename = 'montaged_offsets.txt'
prealigned_offsets_filename = 'prealigned_offsets.txt'
aligned_offsets_filename = 'aligned_offsets.txt'

dir_path = os.path.join(bucket_dir_path, datasets_dir_path, cur_dataset, prealigned_dir_path)
offsets_csv = open(os.path.join(dir_path, prealigned_offsets_filename))
image_reader = csv.reader(offsets_csv, delimiter='\t')

for row in image_reader:
  filename = row[0]
  i_offset = int(row[1])
  j_offset = int(row[2])
  path = os.path.join(dir_path, filename)
  imp = IJ.openImage(path)
  name = imp.getTitle()
  ip = imp.getProcessor()
  ip.setInterpolationMethod(ImageProcessor.BILINEAR)
  ip.translate(j_offset, i_offset)
  ip = ip.resize(int(imp.getWidth()*resize_factor)) # will scale & crop
  new_imp = ImagePlus(name, ip)
  new_imp.show()

# To put to stack:
# Image > Stacks > Images to Stack (align to top-left)