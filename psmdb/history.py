import atexit
import os

historyPath = os.path.expanduser(".pyhistory")

try:
    import readline
except ImportError:
    print("Module readline not available.")
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")
    if os.path.exists(historyPath):
        readline.read_history_file(historyPath)

def save_history(historyPath=historyPath):
    try:
        import readline
    except ImportError:
        print("Module readline not available.")
    else:
        readline.write_history_file(historyPath)
atexit.register(save_history)
