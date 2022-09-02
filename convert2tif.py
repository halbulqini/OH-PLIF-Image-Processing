#import  numpy as np
#import image
import ReadIM
from PIL import Image
import glob
import os
base_dir = r"D:\Hard Disk\Work\_MSC\OH Total Data\IM7\FA-S-07-010-10000_01\Cycle=00003"
new_dir = base_dir +"\\"+ base_dir.split("\\")[-1]
os.mkdir(new_dir)
ls_dir = glob.glob(base_dir + r"\*.im7")
for file_name in ls_dir:
    buffer, atts  = ReadIM.extra.get_Buffer_andAttributeList(file_name)
    v_array, _ = ReadIM.extra.buffer_as_array(buffer)
    im = Image.fromarray(v_array[0])
    new_file_name = new_dir +"\\"+ file_name.split("\\")[-1].replace(".im7", ".tif")
    im.save(new_file_name)
ReadIM.DestroyBuffer(buffer)
ReadIM.DestroyAttributeListSafe(atts)