import sys
import EcoHab

if __name__ == "__main__":
    try:
        path = sys.argv[1]
    except IndexError:
        sys.exit("No path given")

    data = EcoHab.EcoHabData(path)
    
        
