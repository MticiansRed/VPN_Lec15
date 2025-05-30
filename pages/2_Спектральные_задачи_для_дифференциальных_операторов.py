import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Спектральная задача", 
    "Эллиптический оператор",
    "Интегральное тождество",        
    "Вариационная задача", 
    "Конечно-элементное решение",    
    "Матричная задача",                       
    )
)
  
if menu == "Спектральная задача":
    r"""
##### Спектральная задача

В ограниченной области $\Omega$ определим дифференциальный оператор $\mathcal A$

Для достаточно гладких функций $u(x), \ x \in \Omega$, $\\ \quad $ удовлетворяющих однородным граничным условиям на $\partial \Omega$, $\\ \quad$ ищутся нетривиальные решения ($u(x) \not \equiv 0, \ x \in \Omega$) 

$\begin{aligned}
\mathcal A u = \lambda u
\end{aligned}$

Решение спектральной задачи

* $\lambda_k, \quad k = 1,2, \ldots - $ собственные значения
* $\phi_k(x), \quad k = 1,2, \ldots - $ собственные функции

Более общие спектральные задачи

$\begin{aligned}
\mathcal A u = \lambda \mathcal B u
\end{aligned}$

$\mathcal A ,  \mathcal B - $ дифференциальные операторы

    """

if menu == "Эллиптический оператор":
    r"""
##### Эллиптический оператор

Оператор второго порядка

$\begin{aligned} 
 \mathcal A u = - \operatorname{div} \big (k(x) \operatorname{grad} u + \bm a(x) u \big ) +
 \bm b(x) \cdot \operatorname{grad} u + c(x) u ,
 \quad x \in  \Omega 
\end{aligned}$

на множестве функций

$\begin{aligned} 
 k(x) \frac{\partial u}{\partial n} + \sigma (x) u = 0,
 \quad x \in   \partial \Omega
\end{aligned}$

Коэффициенты $k(x) > 0, \ c(x) \geq 0, \ \sigma(x) \geq 0$

Конвективный перенос

* $\bm a(x) ~-~$ консервативный (дивергентный)
* $\bm b(x) ~-~$ характеристический (недивергентный)
    """
    
if menu == "Интегральное тождество":
    r"""
##### Интегральное тождество

$V$ - пространство достаточно гладких функций

Для $u,v \in V$

$\begin{aligned}
(\mathcal{A} u - \lambda u, v) = 0 
\end{aligned}$  

С учетом граничных условий

$\begin{aligned}
& \int_{\Omega} k(x) \operatorname{grad} u \operatorname{grad} v \, d x + \int_{\Omega}  c(x) u \, v \, dx +
\int_{\partial \Omega}  \sigma(x) u \, v \, dx \\ 
& - \int_{\Omega} \operatorname{div} (\bm a(x) u) v \, dx +
\int_{\Omega}  \bm b(x) \cdot \operatorname{grad} u \, v \, d x \\ 
& = \lambda \int_{\Omega} u \, v \, d x 
\end{aligned}$

    """

if menu == "Вариационная задача":
    r"""
##### Вариационная задача

Найти $\lambda, u \in V$ такие, что 

$\begin{aligned}
a(u,v) = \lambda b(u,v),
\quad \forall v \in V
\end{aligned}$ 

Билинейные формы

$\begin{aligned}
a(u,v) & = \int_{\Omega} k(x) \operatorname{grad} u \operatorname{grad} v \, d x + \int_{\Omega}  c(x) u \, v \, dx +
\int_{\partial \Omega}  \sigma(x) u \, v \, dx \\ 
& - \int_{\Omega} \operatorname{div} (\bm a(x) u) v \, dx +
\int_{\Omega}  \bm b(x) \cdot \operatorname{grad} u \, v \, d x 
\end{aligned}$ 

$\begin{aligned}
b(u,v) = \int_{\Omega} u \, v \, d x 
\end{aligned}$ 

    """
if menu == "Конечно-элементное решение":
    r"""
##### Конечно-элементное решение

Конечно-элементный базис в $V_h \subset V$

$\begin{aligned}
 \varphi_i(x),
 \quad i = 1,2, \ldots, n 
\end{aligned}$

$\begin{aligned}
 \varphi_i(x_j) = \left \{ \begin{array}{cc}
  1 ,  &  i = j \\
  0 ,  &  i \neq j \\
\end{array}
\right .
\end{aligned}$

Приближенное решение

$\begin{aligned}
   y(x) = \sum_{i=1}^{n} y_i \varphi_i(x) ,
   \quad y_i = y(x_i),
   \quad i = 1,2, \ldots, n
\end{aligned}$

    """

if menu == "Матричная задача":
    r"""
##### Матричная задача

Вариационная постановка

$\begin{aligned}
y \in V_h = \, ?
\end{aligned}$ 

$\begin{aligned}
a(y,v) = \lambda b(y,v),
\quad \forall v \in V_h
\end{aligned}$ 

Метод Галеркина

$\begin{aligned}
 \sum_{j=1}^{n} a(\varphi_j, \varphi_i)  y_j  = \lambda
 \sum_{j=1}^{n} b(\varphi_j, \varphi_i)  y_j
 \quad i = 1,2, \ldots , n
\end{aligned}$

Спектральная задача линейной алгебры

$\begin{aligned}
A y = \lambda B y
\end{aligned}$ 

Матрицы

$\begin{aligned}
A = \{a_{ij} \}, \quad B = \{b_{ij} \}
\end{aligned}$ 

с элементами

$\begin{aligned}
  a_{ij} = a(\varphi_j, \varphi_i) ,
  \quad b_{ij} = b(\varphi_j, \varphi_i) ,
  \quad i, j = 1,2, \ldots , n
\end{aligned}$


    """



