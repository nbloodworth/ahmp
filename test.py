from modules.mhcdb import MHCdatabase
import sys

loc=sys.argv[1]
db=MHCdatabase(location=loc,allele="A*02:01",toprint="all")