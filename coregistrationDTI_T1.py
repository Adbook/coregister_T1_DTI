import os
import os.path
import shutil
from shutil import copyfile

dtiprocess_path = '/tools/bin_linux64/dtiprocess'
dtiestim_path = '/tools/bin_linux64/dtiestim'
BRAINSFit_path = '/tools/bin_linux64/BRAINSFit'
ANTS_path = '/tools/bin_linux64/ANTS'
ITKTransformTools_path = '/tools/bin_linux64/ITKTransformTools'
ResampleScalarVectorDWIVolume_path = '/tools/bin_linux64/ResampleScalarVectorDWIVolume'
polydatatransform_path = '/tools/bin_linux64/polydatatransform'
polydatamerge_path = '/tools/bin_linux64/polydatamerge'



#structure used to simplify all the differents paths handling
class Case:
  def __init__(self, name = '', path = '', DTIpath = '', DWIpath = '', brainMaskPath = '', innerSurfacesPaths = '', T1path = ''):
    self.name = name
    self.path = path
    self.DTIpath = DTIpath
    self.DWIpath = DWIpath
    self.brainMaskPath = brainMaskPath
    self.innerSurfacesPaths = innerSurfacesPaths
    self.T1path = T1path



  def printself(self):
    print("Name : " + self.name)
    print("Path : " + self.path)
    print("DTI : " + self.DTIpath)
    print("DWI : " + self.DWIpath)
    print("Mask : " + self.brainMaskPath)
    print("lh surf : " + self.innerSurfacesPaths[0][0])
    print("rh surf : " + self.innerSurfacesPaths[0][1])
    print("T1 : " + self.T1path)

  name = ''
  path = ''
  
  #input files
  DTIpath = ''
  DWIpath = ''
  brainMaskPath = ''
  #pairs of surfaces, with lh first and rh second
  innerSurfacesPaths = []
  
  T1path = ''

  #temporary files
  B0path = ''
  ADpath = ''
  FApath = ''

  rigidPath = ''
  affinePath = ''
  affineITKPath = ''
  affineInvITKPath = ''

  regPath = ''

  regWarpPath = ''
  regInvWarpPath = ''

  globalWarpPath = ''
  globalInvWarpPath = ''
  #output
  regT1Path = ''
  regInnerSurfaces = []
  outputMergedInnerSurfacePath = ''


def check_for_file(path, name):
  if os.path.isfile(path):
    print("Found " + name + " " + path)
    return True
  else:
    print("Error: " + name + " not found in " + path)
    return False

def check_case_for_surfaces(case):
  for surf in case.innerSurfacesPaths:
    if not check_for_file(surf[0], surf[0]) or not check_for_file(surf[1], surf[1]):
      return False
  return True

def check_case_for_input_files(case):
  counter = 0
  counter += check_for_file(case.DTIpath, "DTI file")
  counter += check_for_file(case.DWIpath, "DTI file")
  counter += check_for_file(case.brainMaskPath, "Brain mask file")
  counter += check_for_file(case.T1path, "T1 file")
  if not check_case_for_surfaces(case):
    return False
  return counter == 4

def run_command(command):
  print("Running command : " + command)

  if os.system(command) != 0:
    print("Errors were encountered while running this command: \n" + command)
    return False
  return True  


def coregister_T1_to_DTI(case, general_output_path):
  print("Checking data for case " + case.name + "...")
  if check_case_for_input_files(case):
    print("All files were found for case " + case.name)
  else:
    print("Some files were missing for case" + case.name + "... skipping")
    return False

  print("Copying files...") ## we copy all files to a new directory and create a new case that will use the copied files
  output_path = general_output_path + "/" + case.name
  out_dti = output_path + "/" + "DTI"
  out_smri = output_path + "/" + "sMRI"
  if not os.path.isdir(output_path):
    os.mkdir(output_path)
    os.mkdir(out_dti)
    os.mkdir(out_smri)
  
  new_surfaces = []
  for surfaces in case.innerSurfacesPaths:
    new_surf = [out_smri + "/" + case.name + "_" + os.path.basename(surfaces[0]) , out_smri + "/" + case.name + "_" + os.path.basename(surfaces[1])]
    new_surfaces.append(new_surf)

  new_case = Case(case.name, output_path, out_dti + "/" + case.name + "_DTI.nrrd", out_dti + "/" + case.name + "_DWI.nrrd", out_dti + "/" + case.name + "_brainmask.nrrd", new_surfaces, out_smri + "/" + case.name + "_T1.nrrd")

  shutil.copy(case.DTIpath, new_case.DTIpath) #copying files
  shutil.copy(case.DWIpath, new_case.DWIpath)
  shutil.copy(case.brainMaskPath, new_case.brainMaskPath)
  for i in range(0, len(case.innerSurfacesPaths)):  #copying surfaces
    print( case.innerSurfacesPaths[i][0] + "     " + new_case.innerSurfacesPaths[i][1] )
    shutil.copy(case.innerSurfacesPaths[i][0], new_case.innerSurfacesPaths[i][0])
    shutil.copy(case.innerSurfacesPaths[i][1], new_case.innerSurfacesPaths[i][1])
  shutil.copy(case.T1path, new_case.T1path)

  case = new_case

  print("Computing B0 and AD...")
  case.B0path = out_dti + "/" + case.name + "_B0.nrrd"
  case.ADpath = out_dti + "/" + case.name + "_AD.nrrd"
  case.FApath = out_dti + "/" + case.name + "_FA.nrrd"

  command =  dtiestim_path
  command += " --dwi_image " + case.DWIpath
  command += " --B0 " + case.B0path
  command += " --tensor_output to_be_removed.nrrd"
  command += " -M " + case.brainMaskPath
  command += " --correction nearest -m wls -t 0 "
  
  if not run_command(command) or not check_for_file(case.B0path, " B0 "):
    return False

  os.system("rm to_be_removed.nrrd")

  command =  dtiprocess_path
  command += " --dti_image " + case.DTIpath
  command += " --lambda1_output " + case.ADpath
  command += " --fa_output " + case.FApath

  if not run_command(command) or not check_for_file(case.ADpath, " AD ") or not check_for_file(case.FApath, " FA "):
    return False

  print("Registering...")
  case.rigidPath = out_smri + "/" + case.name + "_reg_rigid.txt"
  case.affinePath = out_smri + "/" + case.name + "_regAffine.txt"
  case.affineITKPath = out_smri + "/" + case.name + "_reg_affine_itk.txt"
  case.affineInvITKPath = out_smri + "/" + case.name + "_reg_affine_inv_itk.txt"
  case.regPath = out_smri + "/" + case.name + "_reg.nrrd"
  case.regWarpPath = out_smri + "/" + case.name + "_regWarp.nrrd"
  case.regInvWarpPath = out_smri + "/" + case.name + "_regInverseWarp.nrrd"
  case.regT1Path = out_dti + "/" + case.name + "_T1_regDTI.nrrd"
  case.globalWarpPath = out_smri + "/" + case.name + "_GlobalWarp.nrrd"
  case.globalInvWarpPath = out_smri + "/" + case.name + "_GlobalInvWarp.nrrd"
  case.outputMergedInnerSurfacePath = out_dti + "/" + case.name + "_combined_innerSurf.vtk"

  command = BRAINSFit_path
  command += " --movingVolume " + case.T1path
  command += " --fixedVolume " + case.ADpath
  command += " --linearTransform " + case.rigidPath
  command += " --useRigid --initializeTransformMode useCenterOfHeadAlign"

  if not run_command(command) or not check_for_file(case.ADpath, " AD "):
    return False

  command = ANTS_path + " 3 -r Gauss\[3,1\] -i 100x30x5 -t SyN\[0.25\] "
  command += " -m CC\[" + case.FApath + "," + case.T1path + ",1,8\]"
  command += " -m CC\[" + case.B0path + "," + case.T1path + ",2,8\]"
  command += " -o " + case.regPath 
  command += " --initial-affine " + case.rigidPath
  command += " --use-all-metrics-for-convergence"

  if not run_command(command) or not check_for_file(case.regPath, " reg "):
    return False

  command = ITKTransformTools_path + " MO2Aff " + case.affinePath + " " + case.affineITKPath
  if not run_command(command) or not check_for_file(case.affineITKPath, " affine transform "): 
    return False

  command = ITKTransformTools_path + " invert " + case.affineITKPath + " " +  case.affineInvITKPath
  if not run_command(command) or not check_for_file(case.affineInvITKPath, " Inverse affine transform "): 
    return False

  command = ITKTransformTools_path + " concatenate " + case.globalWarpPath 
  command += " -r " + case.ADpath + " " +  case.affineITKPath + " " +  case.regWarpPath + " displacement "
  if not run_command(command) or not check_for_file(case.globalWarpPath, "global disp field"):
    return False

  command = ITKTransformTools_path + " concatenate " + case.globalInvWarpPath
  command += " -r " + case.T1path + " "  + case.regInvWarpPath + " displacement " +  case.affineInvITKPath  #+ case.affineInvITKPath
  if not run_command(command) or not check_for_file(case.globalInvWarpPath, "global inverse disp field"):
    return False

  command = ResampleScalarVectorDWIVolume_path + " "
  command += case.T1path + " "
  command += case.regT1Path
  command += " -H " + case.globalWarpPath + " --hfieldtype displacement "
  command += " -R " + case.ADpath
  if not run_command(command) or not check_for_file(case.regT1Path, "registered T1 image"):
    return False

  lh_list = ""
  rh_list = ""
  os.system("rm " + case.outputMergedInnerSurfacePath)
  for surfaces in case.innerSurfacesPaths:
    regInnerSurfaces = [out_dti + "/" + case.name + "_" + os.path.basename(surfaces[0]) , out_dti + "/" + case.name + "_" + os.path.basename(surfaces[1])]
    case.regInnerSurfaces.append(regInnerSurfaces) #might not be needed
    lh_list += regInnerSurfaces[0] + " "
    rh_list += regInnerSurfaces[1] + " "
    os.system("rm " + regInnerSurfaces[0])
    os.system("rm " + regInnerSurfaces[1])
    print("Transforming surfaces...")
    command = polydatatransform_path
    command += " --fiber_file " + surfaces[0] 
    command += " -D " + case.globalInvWarpPath + " --inverty --invertx "
    command += " -o " + regInnerSurfaces[0]
    if not run_command(command) or not check_for_file(regInnerSurfaces[0], " transformed surface "):
      return False
    command = polydatatransform_path
    command += " --fiber_file " + surfaces[1] 
    command += " -D " + case.globalInvWarpPath + " --inverty --invertx "
    command += " -o " + regInnerSurfaces[1] 
    if not run_command(command) or not check_for_file(regInnerSurfaces[1], " transformed surface " ):
      return False
    break

  print("Merging...")
  command = polydatamerge_path
  command += " -f " + lh_list
  command += " -g " + rh_list
  command += " -o " + case.outputMergedInnerSurfacePath
  if not run_command(command) or not check_for_file(case.outputMergedInnerSurfacePath, "output registered inner surface"):
    return False

  if not check_for_file(case.outputMergedInnerSurfacePath, "final output"):
    print("Process finished, but the output is not present. Something went wrong.")
    return False
  print("Done !")

  return True



def generate_ADNI_case(path):
  case = Case()
  prefix = 'NO DTI FOUND'
  f = []
  for (dirpath, dirnames, filenames) in os.walk(path + "/DTI"):
    f.extend(filenames)
    break
  for filename in f:
    if filename[0] == 'S' and filename.endswith('.nrrd'):
      prefix = filename.replace('.nrrd', '')
      break

  case.name = prefix
  case.path = path
  case.DTIpath = path + "/DTI/resampling/" + prefix + "_resampled_DTI.nrrd"
  case.DWIpath = path + "/DTI/resampling/" + prefix + "_resampled_dwi.nrrd"
  case.brainMaskPath = path + "/DTI/resampling/" + prefix + "_resampled_brain_mask.nrrd"
  case.innerSurfacesPaths =  [[path + "/vtk/lh.white_labeled.vtk", path + "/vtk/rh.white_labeled.vtk"]]
  case.T1path = path + "/mri/brain.nrrd"
  return case

def generate_all_ADNI_cases(path): #in: path to the FS_Data dir, out: list of cases
  cases = []
  #get list of all dirs
  d = []
  for (dirpath, dirnames, filenames) in os.walk(path):
    d.extend(dirnames)
    break
  #generate cases for every directory (checks will be done later)
  for dirname in d:
    if dirname.startswith('ADNI_'):
      cases.append(generate_ADNI_case(path + '/' + dirname))
  return cases

def main(args):

  #test_case = Case("S177586", "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553", "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/DTI/resampling/S177586_resampled_DTI.nrrd", "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/DTI/resampling/S177586_resampled_dwi.nrrd", "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/DTI/resampling/S177586_resampled_brain_mask.nrrd", [["/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/vtk/lh.white_labeled.vtk", "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/vtk/rh.white_labeled.vtk"]], "/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553/mri/brain.nrrd")
  #test_case_2 = generate_ADNI_case("/ASD/Martin_Data/ADNI/FS_Data/ADNI_005_S_4910_20121213135151553")
  
  
  results_dir = "/ASD/Martin_Data/ADNI/processing/DTI/test"

  #get a list of all the files
  cases = generate_all_ADNI_cases('/ASD/Martin_Data/ADNI/FS_Data')

  #check that all expected files are present, and only retain the valid cases
  valid_cases = []
  for case in cases:
    if check_case_for_input_files(case):
      valid_cases.append(case)
  print("Found " + str(len(valid_cases)) + " valid cases")

  #run the process for every case
  for case in valid_cases:
    coregister_T1_to_DTI(case, results_dir)

if __name__ == "__main__":
    main(0)