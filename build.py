#!/bin/env python3

import sys, os, shutil, stat

def rmminusrf(top):
    if not os.path.isdir(top):
        return
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(top)

name = sys.argv[1]
base = "build/"+name
base1 = sys.argv[2]
files = list(sys.argv[3:])

rmminusrf(base)

from cx_Freeze import setup, Executable
                                                                                                                             
# Dependencies are automatically detected, but it might need fine tuning.
# build_exe_options = {"packages": ["os","encodings.ascii","encodings.utf_8","encodings.utf_16_le","encodings.latin_1","encodings.string_escape","encodings.hex_codec","encodings","_strptime"], "includes": ["xlwt.ExcelFormulaParser"], "excludes": ["tkinter"]}
build_exe_options = {"packages": ["encodings","_strptime"], 
                     "includes": ["xlwt.ExcelFormulaParser",
                                  "xlwt.ExcelFormulaLexer",
                                  "sqlobject.boundattributes",
                                  "sqlobject.declarative",
                                  "secrets",
                                 ], 
                     "excludes": ["tkinter"]}

sys.argv[1:] = ['build']
setup(executables = [Executable(f, base=None) for f in files],
      options     = {"build_exe": build_exe_options })
shutil.move('build/exe.linux-x86_64-3.6',base)

datadir = "%s.data"%(name,)
if os.path.isdir(datadir):
    for f in os.listdir(datadir):
        if os.path.isdir(os.path.join(datadir,f)):
            if f != '.svn':
                shutil.copytree(os.path.join(datadir,f),os.path.join(base,f))
        else:
            shutil.copyfile(os.path.join(datadir,f),os.path.join(base,f))
            if f.endswith('.exe') or f.endswith('.sh'):
                os.chmod(os.path.join(base,f),stat.S_IRWXU|stat.S_IRGRP|stat.S_IXGRP|stat.S_IROTH|stat.S_IXOTH)
            elif f.startswith('lib') and '.so.' in f:
                os.chmod(os.path.join(base,f),stat.S_IRWXU|stat.S_IRGRP|stat.S_IXGRP|stat.S_IROTH|stat.S_IXOTH)
                shutil.move(os.path.join(base,f),os.path.join(base+'/lib',f))
            else:
                os.chmod(os.path.join(base,f),stat.S_IREAD|stat.S_IRGRP|stat.S_IROTH)

for root, dirs, files in os.walk(base):
    dirs1 = []
    for d in dirs:
        if d == '.svn':
            rmminusrf(os.path.join(root,d))
            # os.rmdir(os.path.join(root,d))
        else:
           dirs1.append(d)
    dirs = dirs1

os.system("cd %s; tar -czf %s.tgz %s"%("build",base1,name))
