#!/bin/env python

import sys
import uuid
import numpy as np
import matplotlib.pyplot as plt


path = 'docs/images/' + uuid.uuid4().hex + ".svg"

bw = 0.15
steps = (sys.argv[1:4])
data = sys.argv[4:16]


speedup = [  
    1, float(data[0]) / float(data[3]), float(data[0]) / float(data[6]), float(data[0]) / float(data[9]),
    1, float(data[1]) / float(data[4]), float(data[1]) / float(data[7]), float(data[1]) / float(data[10]),
    1, float(data[2]) / float(data[5]), float(data[2]) / float(data[8]), float(data[2]) / float(data[11])
]

timing = [
    [ steps[0], float(data[0]), float(data[3]), float(data[6]), float(data[9])  ],
    [ steps[1], float(data[1]), float(data[4]), float(data[7]), float(data[10]) ],
    [ steps[2], float(data[2]), float(data[5]), float(data[8]), float(data[11]) ]
]


plt.style.use('ggplot')
plt.figure(figsize=(9, 4))

plt.plot(('1', '2', '4', '8'), speedup[0:4],  alpha=0.5, label=steps[0], linewidth=2, color='red')
plt.plot(('1', '2', '4', '8'), speedup[4:8],  alpha=0.5, label=steps[1], linewidth=2, color='green')
plt.plot(('1', '2', '4', '8'), speedup[8:12], alpha=0.5, label=steps[2], linewidth=2, color='blue')

plt.ylabel("Speedup factor")
plt.xlabel("Nodes")

plt.table(cellText=timing, colLabels=('Steps', '1', '2', '4', '8'), rowLoc="right", loc="bottom", bbox=[0, -0.7, 1, 0.4])
plt.subplots_adjust(bottom=0.4)


plt.title('')
plt.legend()

plt.savefig(path)
#plt.show()
print(path)