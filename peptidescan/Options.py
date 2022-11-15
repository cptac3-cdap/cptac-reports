
# from pkg_resources import resource_stream, resource_exists, resource_filename
try:
    from configparser import SafeConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser
import platform, os, os.path, sys

def escape_spaces(s):
    if " " not in s:
        return s
    if s.startswith('"') and s.startswith('"'):
        return s
    return '"%s"'%s

class Config(object):
    def __init__(self,section="PeptideScan"):

        pkgdir = __file__
        while not os.path.isdir(pkgdir):
            pkgdir,dummy = os.path.split(pkgdir)
        # print >>sys.stderr, pkgdir

        platname = platform.uname()
        archdir='-'.join([platname[0],platname[4]]).strip('-')
        iniName = '%s.ini'%archdir

        cfg = SafeConfigParser()
        # if resource_exists(__name__,iniName):
        #    assert(False)
        # cfg.readfp(resource_stream(__name__,iniName),iniName)
        # else:
        cfg.read(os.path.join(pkgdir,iniName))

        assert cfg.has_section(section), "Can't find section %s in config file: %s."%(section,iniName)
        exeBase = cfg.get(section,"exeBase")
        exeExtn = cfg.get(section,"exeExtn")
        exeName = "%s-%s%s"%(exeBase,archdir,exeExtn)

        if cfg.has_option(section,"exePath"):
            exePath = cfg.get(section,"exePath",vars=dict(pkgDir=pkgdir))
        elif cfg.has_option(section,"exeDir"):
            exePath = os.path.join(cfg.get(section,"exeDir",
                                           vars=dict(pkgDir=pkgdir)),exeName)
        else:
            assert(False)
            # exePath = resource_filename(__name__, exeName)

        # print >>sys.stderr, exePath
        assert(os.path.exists(exePath))

        if cfg.has_option(section,"exeOpt"):
            exeOpt  = cfg.get(section,"exeOpt")
        else:
            exeOpt  = ""
        self.exe = exePath
        self.opt = exeOpt

class Options(object):
    def __init__(self,**kw):
        self.opt = {}
        for k,v in kw.items():
            self.set(k,v)

    def set(self,name,value=None):
        self.opt[name] = value

    def unset(self,name):
        del self.opt[name]

    def get(self,name):
        return self.opt.get(name,None)

    def has(self,name):
        return name in self.opt

    def __str__(self):
        l = []
        for k,v in self.opt.items():
            if v != None:
                l.append("-%s %s"%(k,escape_spaces(str(v))))
            else:
                l.append("-%s"%(k,))
        return " ".join(l)
