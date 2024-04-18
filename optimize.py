import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import functions as fun
from smoothing import *

np.savetxt('test.txt', np.array([0, 1]), fmt='%.1e', delimiter=',')
"""
vitadata1_S = np.array([
    [106.0, 209.0, 76.0], [107.0, 206.0, 77.0], [106.625, 206.0, 82.75], [106.375, 205.125, 86.0],
    [105.0, 205.25, 88.625], [103.5, 204.875, 91.75], [104.375, 203.125, 95.25], [105.375, 202.375, 98.625],
    [106.5, 201.5, 103.0], [108.0, 199.125, 107.125], [107.625, 195.75, 112.625], [105.625, 191.25, 119.75],
    [103.5, 187.125, 125.75], [101.0, 185.625, 131.125], [98.0, 186.125, 134.125], [96.375, 186.0, 135.875],
    [97.0, 184.125, 138.625], [100.0, 180.0, 142.375], [105.125, 175.0, 146.625], [110.375, 170.5, 150.75],
    [114.5, 166.25, 155.625], [117.875, 163.5, 158.75], [121.625, 161.25, 161.625], [126.125, 158.875, 164.75],
    [130.5, 157.625, 164.5], [135.0, 156.75, 162.625], [139.75, 156.625, 159.5], [144.375, 157.5, 155.375],
    [149.5, 157.25, 152.5], [155.125, 156.25, 149.625], [161.5, 155.75, 146.125], [169.25, 154.125, 143.375],
    [177.5, 153.25, 140.625], [186.0, 153.625, 137.875], [194.125, 152.75, 135.5], [200.75, 153.875, 132.75],
    [207.25, 157.25, 129.875], [213.0, 159.875, 127.5], [218.5, 163.5, 125.375], [225.0, 166.75, 122.625],
    [232.375, 168.125, 119.375], [242.0, 169.375, 116.5], [252.375, 170.25, 115.0], [261.0, 172.0, 113.625],
    [267.75, 175.0, 112.0], [274.25, 178.75, 110.625], [280.625, 184.25, 108.25], [286.0, 190.625, 105.875],
    [290.875, 197.0, 104.625], [291.0, 203.75, 105.0], [287.25, 210.25, 106.375], [285.125, 215.75, 107.875],
    [281.875, 221.875, 108.75], [278.5, 229.25, 108.375], [276.5, 236.375, 107.375], [274.625, 242.875, 106.75],
    [272.0, 249.875, 106.25], [268.125, 257.125, 105.625], [264.375, 264.5, 106.0], [259.25, 273.5, 105.875],
    [254.125, 281.5, 107.0], [250.125, 288.5, 109.375], [245.75, 295.75, 110.5], [240.375, 303.125, 111.875],
    [235.625, 311.625, 112.5], [231.125, 320.5, 112.5], [226.375, 329.5, 113.0], [224.125, 336.5, 114.125],
    [224.125, 341.5, 115.375], [224.375, 345.75, 116.0], [222.125, 349.75, 116.875], [218.25, 354.0, 118.5],
    [213.875, 357.125, 122.0], [210.375, 359.875, 126.25], [208.625, 363.0, 129.75], [207.0, 365.125, 134.5],
    [204.875, 367.625, 137.875], [203.125, 370.625, 139.625], [203.0, 372.375, 141.5], [202.5, 374.375, 141.75],
    [201.75, 376.375, 143.125], [201.75, 376.375, 146.375], [201.625, 376.5, 149.0], [201.25, 377.625, 151.875],
    [202.25, 377.5, 155.625], [205.125, 375.0, 160.25], [208.875, 370.625, 166.0], [213.0, 367.0, 169.875],
    [216.25, 364.125, 171.625], [215.875, 363.875, 173.125], [214.375, 364.625, 175.875], [216.875, 363.375, 179.75],
    [221.625, 363.125, 182.875], [228.5, 364.125, 182.625], [237.375, 365.0, 178.875], [245.375, 364.375, 175.75],
    [251.0, 363.875, 172.875], [256.25, 362.625, 172.625], [261.875, 359.625, 177.5], [267.375, 356.125, 183.375],
    [275.0, 351.375, 188.5], [283.875, 347.125, 193.5], [295.5, 343.375, 195.125], [307.875, 342.25, 192.75],
    [318.5, 342.5, 190.25], [329.25, 342.5, 185.625], [337.5, 345.125, 180.0], [344.0, 347.75, 176.875],
    [352.0, 352.125, 173.125], [362.5, 357.25, 170.625], [374.625, 361.0, 170.0], [386.125, 366.75, 167.75],
    [396.75, 370.875, 167.0], [404.25, 373.0, 168.0], [409.0, 376.125, 167.125], [417.0, 377.875, 166.0],
    [427.25, 379.0, 165.25], [437.625, 382.75, 161.875], [447.625, 387.625, 158.125], [454.625, 392.0, 157.375],
    [460.0, 395.375, 158.625], [465.125, 398.5, 160.875], [468.625, 404.75, 163.125], [470.5, 412.375, 164.75],
    [474.125, 419.75, 165.0], [478.875, 427.5, 164.25], [480.5, 435.25, 164.75], [481.5, 443.75, 164.5],
    [481.875, 455.375, 161.25], [480.375, 467.5, 158.625], [481.375, 477.75, 157.875], [482.5, 489.125, 157.25],
    [478.5, 501.125, 158.375], [474.75, 513.25, 158.25], [472.875, 524.375, 156.125], [470.875, 532.25, 156.75],
    [470.125, 539.5, 158.5], [465.625, 548.125, 163.125], [457.5, 557.25, 170.625], [449.25, 568.125, 176.625],
    [440.75, 579.75, 183.0], [433.625, 590.125, 188.0], [429.875, 597.5, 191.875], [428.5, 603.625, 194.5],
    [430.125, 606.75, 197.0], [432.25, 608.625, 201.625], [430.5, 616.625, 203.0], [425.625, 626.5, 204.25],
    [419.375, 636.875, 206.75], [413.375, 647.0, 210.125], [408.0, 652.875, 218.375], [403.25, 658.25, 227.375],
    [398.625, 664.75, 234.75], [395.25, 669.25, 240.5], [394.25, 672.875, 244.625], [394.25, 676.0, 250.0],
    [393.5, 677.125, 258.375], [390.125, 678.5, 267.875], [386.875, 678.25, 277.0], [387.375, 675.375, 285.125],
    [390.875, 673.625, 289.625], [397.25, 672.5, 293.375], [403.375, 672.25, 297.75], [408.125, 672.125, 302.375],
    [414.875, 669.875, 307.5], [422.875, 666.375, 312.75], [432.875, 664.5, 315.25], [442.125, 665.25, 316.125],
    [447.375, 667.25, 317.5], [454.25, 668.375, 317.375], [463.125, 668.125, 319.125], [474.625, 665.375, 323.25],
    [490.0, 662.25, 326.75], [504.625, 662.5, 326.25], [519.625, 664.125, 322.875], [535.25, 665.5, 316.5],
    [550.25, 665.25, 309.5], [567.75, 664.125, 306.75], [585.625, 658.625, 307.875], [604.75, 650.5, 312.125],
    [622.875, 645.75, 314.75], [638.375, 639.375, 316.75], [656.25, 632.875, 321.75], [672.75, 629.5, 327.25],
    [688.125, 626.125, 335.875], [704.0, 625.25, 344.0], [718.875, 626.625, 345.0], [732.5, 628.625, 343.0],
    [746.375, 630.625, 339.125], [761.25, 630.75, 337.0], [771.0, 633.0, 332.0], [791.0, 630.0, 346.0]]).T

vitadata1_S = vitadata1_S[:, :40]
"""
vita_R = s_wm983c_R[:40]
vita_Y = s_wm983c_Y[:40]
vita_G = s_wm983c_G[:40]


def fn_min(params):
    t, r, y, g, m = fun.cc_3_logdelay(params)
    simu_R, simu_Y, simu_G = m + r, y, g
    return np.sum((simu_R - vita_R) ** 2) + np.sum((simu_Y - vita_Y) ** 2) + np.sum((simu_G - vita_G) ** 2)  # *


parameters = np.array([5.101089691028263e+01,
                       4.134421316370464e+00,
                       1.460202315290193e+00,
                       1.253527743481907e+01,
                       7.430237421109924e-01,
                       8.418714315835516e+00,
                       6.563746706018511e-01,
                       8.947126588628151e+00,
                       7.010408140809474e-01,
                       1.163077742270154e+01,
                       4.571958301088177e-01,
                       5.861721344080051e+00,
                       9.650462027760387e-01])

res = minimize(fn_min, parameters, method='Nelder-Mead',
               options={'maxiter': 250})

# bounds=[[0.1,51],[3,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2]]

print(res)
print(res.x)

np.savetxt('results.txt', np.array(res.x), fmt='%.15e', delimiter=',')