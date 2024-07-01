import json, sys

with open (sys.argv[1], 'r') as fh:
    results = json.load (fh)
    ES = 0
    EN = 0
    for b, bl in results["branch attributes"]["0"].items():
        ES += bl ["dS"]
        EN += bl ["dN"]
    
    print ("dS = %g, dN = %g" % ((ES, EN)))
