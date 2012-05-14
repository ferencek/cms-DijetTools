#!/usr/bin/env python

# CSVL
eff0_bb = [0.075190481730750422, 0.092189096135135151, 0.17201223401970492, 0.32047805349296232, 0.44234154305258216]
eff2_bb = [0.52785268605907942, 0.48045368192262933, 0.334041805776981, 0.18690492564265077, 0.11826979776906879]

eff0_qq = [0.78351594923269763, 0.77540241038240232, 0.75544065496584467, 0.7016900412398327, 0.61564715158264538]
eff2_qq = [0.014503055581162171, 0.012011188109140506, 0.016476316213718391, 0.027423668771216252, 0.047029415910807629]

eff0_gg = [0.81152774917400006, 0.80609292179740544, 0.77675438639063388, 0.73871985335251, 0.63226053017346873]
eff2_gg = [0.010221774831434365, 0.010949430552644784, 0.013627544063172525, 0.017896396151228636, 0.039093063936816198]


eff0_sys_up_bb   = [0.010219706183534143, 0.012473030851601145, 0.019620181633986555, 0.040560602922774504, 0.039263984277089076]
eff0_sys_down_bb = [0.0095496487376482506, 0.011704418294451379, 0.018620720162180145, 0.038155264534711394, 0.037553005508526645]
eff2_sys_up_bb   = [0.026973283228365352, 0.027583170855405945, 0.026282110620178281, 0.030936028907094185, 0.020878956880156428]
eff2_sys_down_bb = [0.02630371437099753, 0.026815082188727853, 0.025282860042449895, 0.028531149840199443, 0.019169346359697331]

eff0_sys_up_qq   = [0.016049975148771375, 0.016768161312403978, 0.017323309649634834, 0.03694897121050017, 0.047168475075729799]
eff0_sys_down_qq = [0.015865677350171799, 0.016636813486595087, 0.017143291052891078, 0.035967019734933561, 0.045526424763214904]
eff2_sys_up_qq   = [0.0024042933081408649, 0.0018516552201274818, 0.002533228017061951, 0.0080002323366813576, 0.013864926338439257]
eff2_sys_down_qq = [0.0022212674785785204, 0.0017213973608782146, 0.0023524374455006505, 0.006996947165910916, 0.012191558894962972]

eff0_sys_up_gg   = [0.012108115257458326, 0.011937777396291425, 0.013236654339532921, 0.023335807980699008, 0.038764333554814498]
eff0_sys_down_gg = [0.012017178023226617, 0.011844137863560424, 0.013143485255856792, 0.023047144947170559, 0.037989469621442271]
eff2_sys_up_gg   = [0.0014054360533583413, 0.0014733882028787556, 0.0016892686419348554, 0.0034793318110904799, 0.0093838297543765402]
eff2_sys_down_gg = [0.0013156362179007691, 0.0013806344142037166, 0.0015966813178558063, 0.0031845801147715115, 0.0086009113584815736]

eff0_sys_max_bb = []
eff0_sys_max_rel_bb = []
eff2_sys_max_bb = []
eff2_sys_max_rel_bb = []

eff0_light = []
eff2_light = []
eff0_sys_max_light = []
eff0_sys_max_rel_light = []
eff2_sys_max_light = []
eff2_sys_max_rel_light = []

for i in range(0,len(eff0_bb)):
  eff0_sys_max_bb.append( max( eff0_sys_up_bb[i], eff0_sys_down_bb[i] ) )
  eff0_sys_max_rel_bb.append( max( eff0_sys_up_bb[i], eff0_sys_down_bb[i] ) / eff0_bb[i] )
  eff2_sys_max_bb.append( max( eff2_sys_up_bb[i], eff2_sys_down_bb[i] ) )
  eff2_sys_max_rel_bb.append( max( eff2_sys_up_bb[i], eff2_sys_down_bb[i] ) / eff2_bb[i] )

  eff0_l = 0.5*(eff0_qq[i] + eff0_gg[i])
  eff0_light.append( eff0_l )
  eff0_sys_max_light.append( max( abs(eff0_qq[i] + eff0_sys_up_qq[i] - eff0_l), abs(eff0_qq[i] - eff0_sys_down_qq[i] - eff0_l), abs(eff0_gg[i] + eff0_sys_up_gg[i] - eff0_l), abs(eff0_gg[i] - eff0_sys_down_gg[i] - eff0_l) ) )
  eff0_sys_max_rel_light.append( max( abs(eff0_qq[i] + eff0_sys_up_qq[i] - eff0_l), abs(eff0_qq[i] - eff0_sys_down_qq[i] - eff0_l), abs(eff0_gg[i] + eff0_sys_up_gg[i] - eff0_l), abs(eff0_gg[i] - eff0_sys_down_gg[i] - eff0_l) ) / eff0_l )
  eff2_l = 0.5*(eff2_qq[i] + eff2_gg[i])
  eff2_light.append( eff2_l )
  eff2_sys_max_light.append( max( abs(eff2_qq[i] + eff2_sys_up_qq[i] - eff2_l), abs(eff2_qq[i] - eff2_sys_down_qq[i] - eff2_l), abs(eff2_gg[i] + eff2_sys_up_gg[i] - eff2_l), abs(eff2_gg[i] - eff2_sys_down_gg[i] - eff2_l) ) )
  eff2_sys_max_rel_light.append( max( abs(eff2_qq[i] + eff2_sys_up_qq[i] - eff2_l), abs(eff2_qq[i] - eff2_sys_down_qq[i] - eff2_l), abs(eff2_gg[i] + eff2_sys_up_gg[i] - eff2_l), abs(eff2_gg[i] - eff2_sys_down_gg[i] - eff2_l) ) / eff2_l )

print "eff0_sys_max_bb:"
print eff0_sys_max_bb
print "eff0_sys_max_rel_bb:"
print eff0_sys_max_rel_bb
print "eff2_sys_max_bb:"
print eff2_sys_max_bb
print "eff2_sys_max_rel_bb:"
print eff2_sys_max_rel_bb
print ""
print "eff0_light:"
print eff0_light
print "eff2_light:"
print eff2_light
print ""
print "eff0_sys_max_light:"
print eff0_sys_max_light
print "eff0_sys_max_rel_light:"
print eff0_sys_max_rel_light
print "eff2_sys_max_light:"
print eff2_sys_max_light
print "eff2_sys_max_rel_light:"
print eff2_sys_max_rel_light
