import pandas as pd
import numpy as np
import re
import math

#EvaluaciÃÂ³n REGREX
def evaluate_Fx(str_equ, valX):
  x = valX
  exp, sin, cos = math.exp, math.sin, math.cos
  strOut = re.sub(r"e\s*\^", "exp", str_equ)
  strOut = re.sub(r"(?<=\d)x", '*(x)', strOut).replace("^", "**").strip()
  print(strOut)
  out = eval(strOut)
  return out

def evaluate_Fx_y(str_equ, valX, valY):
  x = valX
  y = valY
  strOut = str_equ.replace("x", '*(x)').replace("y", "*(y)").replace("^", "**")
  out = eval(strOut)
  return out


#Deferencias finitas para derivadas progresivas
def evaluate_derivate_fx(str_equ, x, h):
  x = float(x)
  h = float(h)
  strOut = str_equ.replace("x", '*(x + h)')
  strOut = strOut.replace("^", "**")
  strOut = "-4*(" + strOut + ")"
  out = eval(strOut)
  
  strOut = str_equ.replace("x", '*(x + 2*h)')
  strOut = strOut.replace("^", "**")
  out = out + eval(strOut)
  
  strOut = str_equ.replace("x", '*(x)')
  strOut = strOut.replace("^", "**")
  strOut = "3*(" + strOut + ")"
  out = out + eval(strOut)
  
  out = -out/(2*h)
  return out


def f_transform(equ_str, var, x , n): 
  x_h = f"* ({var} + ({n}) * h)"
  return equ_str.replace(var, x_h).replace("^", "**")
  
def forward_difference(equ_str, x, h):
  out = (4 * evaluate_Fx(equ_str, x+h) - evaluate_Fx(equ_str, x+ 2 *h) - 
          3 * evaluate_Fx(equ_str, x)) / (2 * h)
  return out

def central_difference(equ_str, x, h):
  out = (evaluate_Fx(equ_str, x+h) - evaluate_Fx(equ_str, x-h)) / (2*h)
  return out

def central_difference2(equ_str, x, h):
  out = (evaluate_Fx(equ_str, x - 2*h) - 8 * evaluate_Fx(equ_str, x-h) 
            + 8 * evaluate_Fx(equ_str, x+h) - evaluate_Fx(equ_str, x+2*h)) / (12 * h)
  return out

def forward_difference_2(equ_str, x, y, h):
  out_x = (4 * evaluate_Fx_y(equ_str, x+h, y) - evaluate_Fx_y(equ_str, x+2*h, y) - 
          3 * evaluate_Fx_y(equ_str, x, y)) / (2 * h)
  out_y = (4 * evaluate_Fx_y(equ_str, x, y + h) - evaluate_Fx_y(equ_str, x, y +2 *h) - 
          3 * evaluate_Fx_y(equ_str, x, y)) / (2 * h)
  return (out_x, out_y)

def central_difference_2(equ_str, x, y, h):
  out_x = (evaluate_Fx_y(equ_str, x+h, y) - evaluate_Fx_y(equ_str, x-h, y)) / (2*h)
  out_y = (evaluate_Fx_y(equ_str, x, y + h) - evaluate_Fx_y(equ_str, x, y - h)) / (2*h)
  return (out_x, out_y)

def central_difference2_2(equ_str, x, y, h):
  out_x = (evaluate_Fx_y(equ_str, x-2*h, y) - 8 * evaluate_Fx_y(equ_str, x-h, y) 
            + 8 * evaluate_Fx_y(equ_str, x+h, y) - evaluate_Fx_y(equ_str, x+2*h, y)) / (12 * h)
  out_y = (evaluate_Fx_y(equ_str, x, y - 2*h) - 8 * evaluate_Fx_y(equ_str, x, y - h) 
            + 8 * evaluate_Fx_y(equ_str, x, y + h) - evaluate_Fx_y(equ_str, x, y + 2*h)) / (12 * h)
  return (out_x, out_y)


def evaluate_finite_difference(equ_str, x, h):
  x = float(x)
  h = float(h)
  equ_str = equ_str.lower()
  return pd.DataFrame({
              'f':[equ_str], 'X': [x], 
              'H': [h], 
              'DifFinCen': [central_difference(equ_str, x ,h)],
              'DifFinPro': [forward_difference(equ_str, x, h)],
              'DifFinCen2': [central_difference2(equ_str, x, h)]})
              
def evaluate_finite_difference_2(equ_str, x, y, h):
  x = float(x)
  y = float(y)
  h = float(h)
  return pd.DataFrame({
              'f':[equ_str], 'X': [x], 'Y': [y], 'H': [h], 
              'DifFinCen': ["{:.5f},{:.5f}".format(*central_difference_2(equ_str, x, y ,h))],
              'DifFinPro': ["{:.5f},{:.5f}".format(*forward_difference_2(equ_str, x, y, h))],
              'DifFinCen2': ["{:.5f},{:.5f}".format(*central_difference2_2(equ_str, x, y, h))]})
  
  
  

#Resolverdor de Newton
def newtonSolverX(x0, f_x, eps):
  x0 = float(x0)
  eps = float(eps)
  xn = x0
  error = 1
  arrayIters = []
  arrayXn = []
  arrayErr = []
  
  i = 0
  h = 0.000001
  while(error > eps):
    print("...")
    x_n1 = xn - (evaluate_Fx(f_x, xn)/evaluate_derivate_fx(f_x, xn, h))
    error = abs(x_n1 - xn)
    i = i + 1
    xn = x_n1
    arrayIters.append(i)
    arrayXn.append(xn)
    arrayErr.append(error)

  print("Finalizo...")
  TableOut = pd.DataFrame({'Iter':arrayIters, 'Xn':arrayXn, 'Error': arrayErr})
  return TableOut



def newton_raphson_solver(f_x, x0, k_max, eps):
  xn, k_max, eps = float(x0), float(k_max), float(eps)
  k = 0
  error = float("inf")
  h = 0.000001
  error_vals = []
  xn_vals = []
  while k < k_max and error > eps:
    res = evaluate_Fx(f_x, xn)
    error = abs(res)
    xn_vals.append(xn)
    error_vals.append(error)
    k += 1
    xn = xn - (res / central_difference2(f_x, xn, h))
  return pd.DataFrame({'K': range(k), 'Xn': xn_vals, 'Error': error_vals})


def bisection_solver(f_x, a, b, k_max, eps):
  a, b, k_max, eps = float(a), float(b), float(k_max), float(eps)
  xn = (a + b) / 2
  k = 0
  error = float("inf")
  error_vals = []
  xn_vals = []
  while k < k_max and error > eps:
    res = evaluate_Fx(f_x, xn)
    error = abs(res)
    xn_vals.append(xn)
    error_vals.append(error)
    k += 1
    if evaluate_Fx(f_x, a) * res < 0:
      b = xn
    else:
      a = xn
    xn = (a + b) / 2
  return pd.DataFrame({'K': range(k), 'Xn': xn_vals, 'Error': error_vals})


def quadratic_gd_solver(q_str, c_str, x0_str, e_str, n_str, lr_str, lr_method="constante"):
  q = np.array(eval(q_str))
  c = np.array(eval(c_str))
  x0 = np.array(eval(x0_str))
  e = float(e_str)
  n = float(n_str)
  lr = float(lr_str) 
  
  xn_vals = []
  pk_vals = []
  rate_change_vals = []
  k = 0
  x = x0
  gradient = np.matmul(q, x) + c
  rate_of_change = np.linalg.norm(gradient)

  print(lr_method)
  while rate_of_change > e and k < n:
    # x1 = x0 - a * f'(x0)
    k += 1
    if lr_method == 'constante':
      alpha = lr
    elif lr_method == 'variable':
      alpha = 1/k
    else: 
      alpha = (rate_of_change ** 2) / (np.matmul(np.matmul(gradient.T, q), gradient))
    x = x - alpha * gradient
    gradient = np.matmul(q, x) + c
    rate_of_change = np.linalg.norm(gradient)
    pk_vals.append(str(gradient))
    rate_change_vals.append(str(rate_of_change))
    xn_vals.append(str(x))
  return pd.DataFrame({'K': range(1, k + 1), 'Xn': xn_vals, 'Pk': pk_vals, 'MaxRateChange': rate_change_vals})
  

def ronsenbrock_gd_solver(x0_str, e_str, n_str, lr_str):
  def calc_gradient(x):
    g0 = -400 * (x[0] * x[1] - x[0]**3) - 2 *(1 - x[0])
    g1 = 200 * (x[1] - x[0] ** 2)
    return np.array([g0, g1])

  x0 = np.array(eval(x0_str))
  e = float(e_str)
  n = float(n_str)
  lr = float(lr_str)

  xn_vals = []
  pk_vals = []
  rate_change_vals = []
  k = 0
  x = x0
  gradient = calc_gradient(x)
  rate_of_change = np.linalg.norm(gradient)
  while rate_of_change > e and k <= n:
    k += 1
    x = x - lr * gradient
    gradient = calc_gradient(x)
    rate_of_change = np.linalg.norm(gradient)
    pk_vals.append(str(gradient))
    rate_change_vals.append(str(rate_of_change))
    xn_vals.append(str(x))

  return pd.DataFrame({'K': range(k), 'Xn': xn_vals, 'Pk': pk_vals, 'MaxRateChange': rate_change_vals})
