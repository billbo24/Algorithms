# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:45:44 2020

@author: dingl
"""

#mousnow program
import pyautogui


pyautogui.PAUSE = 1
pyautogui.FAILSAFE = True

print('Press Ctrl-C to quit.')
try:
   while True:
       x, y = pyautogui.position()
       positionStr = 'X: ' + str(x).rjust(4) + ' Y: ' + str(y).rjust(4)

except KeyboardInterrupt:
   print('\nDone.')


print(positionStr, end='')
print('\b' * len(positionStr), end='', flush=True)