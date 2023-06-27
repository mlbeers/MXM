import sympy as sym
from sympy import solve
import numpy as np
from scipy.integrate import nquad
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
import os

class Section:
    
    def __init__(self):
        self.points = []
        self.vec = None
        self.dx = 0.005

    def add_point(self, point):
        self.points.append(point)
        
    def line(self, point1, point2):
        slope = (point2[1] - point1[1])/(point2[0]- point1[0])
        intercept = point1[1] - slope*point1[0]
        if slope < 2 * self.dx and slope > -2 * self.dx:
            slope = 0
        if intercept >= 1 - 2*self.dx and intercept <= 1 + 2*self.dx:
            intercept = 1
        return [slope, intercept]
    
    def top(self):
        if self.points[1][0] >= 1 - self.dx:
            return self.line(self.points[0], self.points[1])
        elif len(self.points) == 4:
            return [self.line(self.points[0], self.points[1]), self.line(self.points[1], self.points[2])]
        else:
            raise ValueError("Section is being described by more than 4 points")
    
    def bottom(self):
        if len(self.points) <= 3:
            return self.line(self.points[2], self.points[0])
        if self.points[3][0] >= 1 - self.dx:
            return self.line(self.points[3], self.points[0])
        elif len(self.points) == 4:
            return [self.line(self.points[2], self.points[3]), self.line(self.points[3], self.points[0])]
        else:
            raise ValueError("Section is being described by more than 4 points")
    
    def top_eq(self):
        x = sym.Symbol('x')
        eqs = []
        if type(self.top()[0]) != list:
            return [self.top()[0]*x + self.top()[1]]
        for eq in self.top():
            eqs.append(eq[0]*x + eq[1])
        return eqs
    
    def bottom_eq(self):
        x = sym.Symbol('x')
        eqs = []
        if type(self.bottom()[0]) != list:
            return [self.bottom()[0]*x + self.bottom()[1]]
        for eq in self.bottom():
            eqs.append(eq[0]*x + eq[1])
        return eqs
    
    def t(self):
        x, y, t = sym.symbols('x y t')
        Mab = np.array([[x, y], [0, 1/x]])
        horo = np.array([[1, 0], [-t, 1]])
        a = horo@(Mab@self.vec)
        return solve(a[1][0], t)[0]
    
    def y(self):
        x, y, t = sym.symbols('x y t')
        Mab = np.array([[x, y], [0, 1/x]])
        horo = np.array([[1, 0], [-t, 1]])
        a = horo@(Mab@self.vec)
        return solve(a[1][0], y)[0]
    
    def t_tan(self):
        times = []
        from sympy.abc import x, t
        bottom = self.bottom_eq()
        for item in bottom:
            d_bottom = sym.diff(item)
            y = self.y()
            d_y = sym.diff(y, x)
            equations = [item - y, d_bottom - d_y]
            solutions = solve(equations, x, t, dict=True)
            times.append(solutions[0][t])
        if len(times) != 1 and abs(times[0] - times[1]) <= 2*self.dx:
            return times[0]
        else:
            return times
        
    def t_enter(self):
        x, y = sym.symbols('x y')
        point = 1
        t = self.t()
        if len(self.top_eq()) == 2:
            point = 2
        return t.subs([(x, self.points[point][0]), (y, self.points[point][1])]) 
            
    def t_edge(self):
        if self.points[0][0] == self.dx:
            return 
        x, y = sym.symbols('x y')
        t = self.t()
        return t.subs([(x, self.points[0][0]), (y, self.points[0][1])]) 
    
    def t_point(self):
        if len(self.bottom_eq()) == 1:
            return None
        x, y = sym.symbols('x y')
        t = self.t()
        return t.subs([(x, self.points[3][0]), (y, self.points[3][1])]) 
    
    def times(self):
        times = []
        if abs(self.t_enter() - 1) > 2*self.dx:
            times.append(1)
        
        t = [self.t_enter(), self.t_tan(), self.t_point(), self.t_edge()]
        for item in t:
            if item != None:
                if type(item) != list:
                    times.append(item)
                else:
                    times.extend(item)
        times.append(10)
        return times
    
    def output(self):
        print("function of y: " + str(self.y()))
        print("times: " + str(self.times()))
        print("top: " + str(self.top_eq()))
        print("bottom: " + str(self.bottom_eq()))
        if len(self.bottom_eq()) == 2:
            print("point_x: " + str(self.points[3]))
        