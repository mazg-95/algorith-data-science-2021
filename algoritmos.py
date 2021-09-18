import pandas as pd
import numpy as np
import re
import math
import matplotlib.pyplot as plt
from sklearn import metrics
plt.rcParams["figure.figsize"] = (20,10)
np.random.seed(5)


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


# Lab 4
def generate_and_save_data(d_str, n_str, path=None):
  print("Generate Data")
  d = int(d_str) #cantidad de columnas para el dataset.
  n = int(n_str) #cantidad de observaciones para el dataset.
  A = np.random.normal(0,1, size=(n,d))
  x_true = np.random.normal(0,1, size=(d,1))
  b = A.dot(x_true) + np.random.normal(0,0.5,size=(n,1))
  if not path:
    path = ''
  path += 'data.npy'
  with open(path, 'wb') as f:
    np.save(f, A)
    np.save(f, x_true)
    np.save(f, b)
  return A, x_true, b

def read_data(path):
  with open(path, 'rb') as f:
    A = np.load(f)
    x_true = np.load(f)
    b = np.load(f)
  return A, x_true, b


### Solucion Cerrada

def compute_error(A, x, y):
  return np.power(np.linalg.norm(np.matmul(A, x) - y), 2)

def closed_solution_solver(path):
  A, x_true, b = read_data(path)
  return closed_solution(A, b, x_true)

def closed_solution(A, b, x_true):
  x = np.matmul(np.matmul(np.linalg.inv(np.matmul(A.T, A)), A.T),b)
  res = compute_error(A, x, b)
  error = metrics.mean_squared_error(x_true, x)
  return pd.DataFrame({"X": str(x_true.reshape(-1).tolist()), 
                       "X*": str(x.reshape(-1).tolist()), 
                       "SumSquaresError": [res],
                       "Error": [error]})


## Gradient Descent

def gd_solver(e_str, n_str, lr_str,batch_size_str="", path="", method="gd"):
  A, x_true, b = read_data(path)
  e = float(e_str)
  n = int(n_str)
  lr = float(lr_str)
  if method == "gd":
    return gradient_descent(A, x_true, b, e, n, lr)
  elif method == "sgd":
    return stochastic_gradient_descent(A, x_true, b, e, n, lr)
  elif method == "mbgd":
    batch_size = int(batch_size_str)
    return mb_gradient_descent(A, x_true, b, e, n, lr, batch_size=batch_size)
  else:
    return gradient_descent(A, x_true, b, e, n, lr)

def compute_f(A, x, y):
  return np.power(np.linalg.norm(np.matmul(A, x) - y), 2)

def compute_e(x, x_true):
  return np.linalg.norm(x - x_true)

def gradient_descent(A, x_true, b, e=10**-8, epochs=20, lr=0.00005):
  ### GD
  x = np.zeros((x_true.shape[0],1))
  k = 0
  gradient = np.matmul(A.T, (np.matmul(A, x) - b))
  error = compute_e(0, gradient)
  f = compute_f(A, x, b)
  x_k = [np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6)]
  f_x = [f]
  errors = [error]
  while k < epochs and e < error: 
    x -= lr * gradient
    gradient = np.matmul(A.T, (np.matmul(A, x) - b))
    f = compute_f(A, x, b)
    error = compute_e(0, gradient)
    x_k.append(np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6))
    f_x.append(f)
    errors.append(error)
    k += 1
  print(f, error)
  return pd.DataFrame({"K": range(len(f_x)), "X*": x_k, 
                        "SumSquaresError": f_x,
                        "Error": errors})



## Stochastic Gradient Descent

def stochastic_gradient_descent(A, x_true, b, e=10**-8, epochs=1000, lr=0.00005):
  ### SGD
  x = np.zeros((x_true.shape[0],1))
  gradient = np.matmul(A.T, (np.matmul(A, x) - b))
  error = compute_e(0, gradient)
  f = compute_f(A, x, b)
  x_k = [np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6)]
  f_x = [f]
  errors = [error]
  k = 0
  while k < epochs and e < error: 
    k += 1
    p = np.random.permutation(A.shape[0])
    A = A[p]
    b = b[p]
    for i in range(A.shape[0]):
      example = A[i,:].reshape(-1,1)
      gradient = np.matmul(example, (np.matmul(example.T, x) - b[i]))
      x -= lr * gradient

      ### SAVE PROGRESS
      f = compute_f(A, x, b)
      error = compute_e(0, np.matmul(A.T, (np.matmul(A, x) - b)))
      x_k.append(np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6))
      f_x.append(f)
      errors.append(error)
    
  print(f, error)
  return pd.DataFrame({"K": range(len(f_x)), "X*": x_k, 
                        "SumSquaresError": f_x,
                        "Error": errors})


## Mini Batch Gradient Descent

def mb_gradient_descent(A, x_true, b, e= 10**-8, epochs=50, lr=0.00005, batch_size = 25):
  ### GD
  n = A.shape[0] // batch_size
  x = np.zeros((x_true.shape[0],1))
  gradient = np.matmul(A.T, (np.matmul(A, x) - b))
  error = compute_e(0, gradient)
  f = compute_f(A, x, b)
  x_k = [np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6)]
  f_x = [f]
  errors = [error]
  k = 0
  while k < epochs and e < error: 
    k += 1 
    p = np.random.permutation(A.shape[0])
    A = A[p]
    b = b[p]

    # mb update
    for i in range(n):
      start, end = i * batch_size, (i+1) * batch_size
      a_mb = A[start:end,]
      b_mb = b[start:end,]
      gradient = np.matmul(a_mb.T, (np.matmul(a_mb, x) - b_mb))
      x -= lr * gradient
      ### SAVE PROGRESS
      f = compute_f(A, x, b)
      error = compute_e(0, np.matmul(A.T, (np.matmul(A, x) - b)))
      x_k.append(np.array2string(x.reshape(-1), separator=',', suppress_small=True, threshold=6))
      f_x.append(f)
      errors.append(error)
  print(f, error)
  return pd.DataFrame({"K": range(len(f_x)), "X*": x_k, 
                        "SumSquaresError": f_x,
                        "Error": errors})


# Problema 2

### GD con Backtracking Line Search

def compute_rosenbrock(x):
  return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2

def backtraking_search(x, gradient, pk):
  lr = 1
  c = 10 ** -4
  rho = 0.5
  def check_armijo_condition(x, lr, gradient, pk):
    return compute_rosenbrock(x + lr * pk) <= (compute_rosenbrock(x) + c * lr * np.matmul(gradient, pk))
  while not check_armijo_condition(x, lr, gradient, pk):
    lr *= rho
  return lr

def ronsenbrock_gd_backtracking_solver(x0_str, e_str="1e-8", n_str = "1000", lr_str=None):
  def calc_gradient(x):
    g0 = -400 * (x[0] * x[1] - x[0]**3) - 2 *(1 - x[0])
    g1 = 200 * (x[1] - x[0] ** 2)
    return np.array([g0, g1])

  x0 = np.array(eval(x0_str))
  e = float(e_str)
  n = int(n_str)
  lr = lr_str and float(lr_str)

  xn_vals = []
  pk_vals = []
  rate_change_vals = []
  k = 0
  x = x0
  step_sizes=[]
  gradient = calc_gradient(x)
  rate_of_change = np.linalg.norm(gradient)
  while rate_of_change > e and k <= n:
    k += 1
    alpha = lr if lr else backtraking_search(x, gradient, -1 * gradient)
    
    xn_vals.append(str(x))
    step_sizes.append(alpha)
    pk_vals.append(str(gradient))
    rate_change_vals.append(str(rate_of_change))
    
    x -= alpha * gradient
    gradient = calc_gradient(x)
    rate_of_change = np.linalg.norm(gradient)

  xn_vals.append(str(x))
  step_sizes.append(0)
  pk_vals.append(str(None)) 
  rate_change_vals.append(rate_of_change)
    

  return pd.DataFrame({'K': range(k+1), 'Lr': step_sizes, 'Xn': xn_vals, 'Pk': pk_vals, 'MaxRateChange': rate_change_vals})



### Metodo de Newton con Backtracking 


def ronsenbrock_newton_backtracking_solver(x0_str, e_str="1e-8", n_str = "3000", lr_str=None):
  def calc_gradient(x):
    g0 = -400 * (x[0] * x[1] - x[0]**3) - 2 *(1 - x[0])
    g1 = 200 * (x[1] - x[0] ** 2)
    return np.array([g0, g1])
  def calc_hessian(x):
    return np.array([[-400*(x[1] - 3 * x[0] ** 2) + 2,-400*x[0]],[-400*x[0],200]])

  x0 = np.array(eval(x0_str))
  e = float(e_str)
  n = int(n_str)
  lr = lr_str and float(lr_str)
  
  k = 0
  x = x0
  gradient = calc_gradient(x)
  h = np.linalg.inv(calc_hessian(x))
  rate_of_change = np.linalg.norm(gradient)
  xn_vals = []
  pk_vals = []
  rate_change_vals = []
  step_sizes=[]

  while rate_of_change > e and k <= n:
    k += 1
    pk = -1 * np.matmul(h.T, gradient)
    alpha = lr if lr else backtraking_search(x, gradient, pk)
    
    xn_vals.append(str(x))
    step_sizes.append(alpha)
    pk_vals.append(str(pk))
    rate_change_vals.append(str(rate_of_change))
    
    x += alpha * pk
    gradient = calc_gradient(x)
    h = np.linalg.inv(calc_hessian(x))
    rate_of_change = np.linalg.norm(gradient)
  
  xn_vals.append(str(x))
  step_sizes.append(0)
  pk_vals.append(str(None)) 
  rate_change_vals.append(rate_of_change)

  return pd.DataFrame({'K': range(k+1), 'Lr': step_sizes, 'Xn': xn_vals, 'Pk': pk_vals, 'MaxRateChange': rate_change_vals})
