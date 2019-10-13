# Rob's example of a class:

class Test:
    def __init__(self, width, height):
        self.width=width
        self.height=height
        self.a=None
    def area(self):
        self.a=self.width*self.height
        return(self.width*self.height)

x = Test(5,3)
print(x.width)
print(x.a)
print(x.area())
print(x.a)



# Revision: Re-write the square class:

class Rectangle:
    def __init__(self, width, height):
        self.width = width
        self.height = height
    def area(self):
        self.a = self.width * self.height
        return self.a

trying = Rectangle(6, 8)
print (trying.area())

# ----------- How to raise an error?

ok = 2
if ok == 1:
    raise Exception("MisMatch", "Bases are not matching")
"""
import Alignment

x=mappy.alignment("Filename")

"""
