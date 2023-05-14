import optiPython as oP
import json
import sys


triFan = oP.triangleFan() # initialize a triangle fan

# JSON string
triInfo = '{"x0": [9.01006524, -4.739905], "T0": 25.887855886616833, "grad0": [0.62334049, 0.78195053], "x1": [ 8.91006524, -4.539905  ], "T1": 25.995574287564278, "grad1": [0.4539905 , 0.89100652], "xHat": [ 9.53879533, -3.97683432], "listIndices": [1.0, 1.0, 1.0, 1.452, 1.0], "listxk": [ [ 9.01006524, -4.739905  ], [ 8.91006524, -4.539905  ], [ 9.23879533, -3.82683432], [ 9.53879533, -3.97683432] ], "listB0k": [ [-0.1,  0.2], [0.22873008, 0.91307067], [0.52873008, 0.76307067] ], "listBk": [ [-0.1,  0.2], [0.22873008, 0.91307067], [0.52873008, 0.76307067] ], "listBkBk1": [ [0.4022869, 0.7895325], [0.33910078, 0.8186617 ], [ 0.3 , -0.15], [ 0.3 , -0.15] ], "plotBefore": 1, "plotAfter": 1, "plotOpti": 1 } '

# Read input from C
triInfo = sys.stdin.read()

# Call function with arguments
stringOut = triFan.outputJSON(triInfo)


# serialize output data as JSON and write to stdout
dataOut = json.dumps(stringOut)
sys.stdout.write(dataOut)


