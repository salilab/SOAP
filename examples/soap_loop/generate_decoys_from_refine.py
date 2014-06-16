from SOAP.loop import *



a=sprefine()
a.dir='/bell2/gqdong/runs/'
a.rundir='60099'
a.nrpl=0
a.nors=175
a.assess_method='SOAP'
a.task=1
a.generate_refine_decoys('lr1')
#a.afterprocessing()
#a.analyze_loop_modeling()
