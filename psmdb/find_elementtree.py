
try:
    import xml.etree.cElementTree as ET
except ImportError:
    try:
        import xml.etree.ElementTree as ET
    except:
        import elementtree.ElementTree as ET
