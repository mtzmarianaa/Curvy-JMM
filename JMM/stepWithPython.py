import optiPython as oP
import json
import sys
import matplotlib.pyplot as plt
from pathlib import Path

# # JSON string
# triInfo = '{"x0": [9.01006524, -4.739905], "T0": 25.887855886616833, "grad0": [0.62334049, 0.78195053], "x1": [ 8.91006524, -4.539905  ], "T1": 25.995574287564278, "grad1": [0.4539905 , 0.89100652], "xHat": [ 9.53879533, -3.97683432], "listIndices": [1.0, 1.0, 1.0, 1.452, 1.0], "listxk": [ [ 9.01006524, -4.739905  ], [ 8.91006524, -4.539905  ], [ 9.23879533, -3.82683432], [ 9.53879533, -3.97683432] ], "listB0k": [ [-0.1,  0.2], [0.22873008, 0.91307067], [0.52873008, 0.76307067] ], "listBk": [ [-0.1,  0.2], [0.22873008, 0.91307067], [0.52873008, 0.76307067] ], "listBkBk1": [ [0.4022869, 0.7895325], [0.33910078, 0.8186617 ], [ 0.3 , -0.15], [ 0.3 , -0.15] ], "plotBefore": 1, "plotAfter": 1, "plotOpti": 1 } '


with open(Path(sys.argv[1])) as f:
    triInfo = f.readline()
    print(triInfo)
    params_dict = json.loads(triInfo)
    nRegions = len(params_dict['listxk']) - 2
    triFan = oP.triangleFan(nRegions) # initialize a triangle fan
    dict_out = triFan.outputJSON(triInfo)
    flagLagrange = triFan.setLagrangeFlag()
    # try:
    #     dict_out = triFan.outputJSON(triInfo)
    # except:
    #     import ipdb; ipdb.set_trace()
    print(flagLagrange, ", ", dict_out["THat"], ", ", dict_out["gradHat"][0], ", ", dict_out["gradHat"][1])
    #str_out = json.dumps(dict_out)
    #print(str_out)



# # Read input from C
# triInfo = open('/Users/marianamartinez/Documents/Curvy-JMM/JMM/update.json', 'r')

# # Call function with arguments
# dict_out = triFan.outputJSON(triInfo)


# # serialize output data as JSON and write to stdout
# with open('/Users/marianamartinez/Documents/Curvy-JMM/JMM/updateSolve.json', 'w') as f:
#     json.dump(dict_out, f)

# if(dict_out["plotAfter"] == 1):
#     plt.show()



