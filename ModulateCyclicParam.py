import math

def ModulateCyclicParam (cyclicParam, cyclicity, shift = 0.0):
    if cyclicParam >= (0  + shift):
        return (cyclicParam - shift) % (cyclicity) + shift 
    else:
        return ModulateCyclicParam (cyclicParam + cyclicity, cyclicity) 
