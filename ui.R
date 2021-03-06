library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
dashboardPage(
    dashboardHeader(title = "Algoritmos en la Ciencia de Datos"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Ceros", tabName = "Ceros"),
            menuItem("Derivación Numerica", tabName = "Derivacion"),
            menuItem("Derivacion Numerica 2", tabName = "Derivacion2"),
            menuItem("Metodo Newton", tabName = "newton_method"),
            menuItem("Metodo Biseccion", tabName = "bisection_method"),
            menuItem("GD QP", tabName = "GD_Q"),
            menuItem("GD Rosenbrock", tabName = "GD_Rosenbrock"),
            menuItem("Gen Data", tabName = "gen_data"),
            menuItem("Closed Solution", tabName="Closed_Solution"),
            menuItem("Batch GD", tabName="bgd_tab"),
            menuItem("SGD", tabName="sgd_tab"),
            menuItem("MiniBatch GD", tabName="mbgd_tab"),
            menuItem("GD Backtracking Rosenbrock", tabName = "GD_Rosenbrock2"),
            menuItem("Newton Backtracking Rosenbrock", tabName = "Newton_Rosenbrock")
        )
    ),
    dashboardBody(
        tabItems(
            tabItem("Ceros",
                    h1("Método de Newton"),
                    box(textInput("ecuacion", "Ingrese la Ecuación"),
                        textInput("initVal", "X0"),
                        textInput("Error", "Error")),
                    actionButton("nwtSolver", "Newton Solver"),
                    tableOutput("salidaTabla")),
            
            tabItem("Derivacion",
                    h1("Diferencias Finitas"),
                    box(textInput("difFinEcu", "Ingrese la Ecuación"),
                    textInput("valorX", "x"),
                    textInput("valorH", "h")),
                    actionButton("diferFinEval", "Evaluar Derivada"),
                    tableOutput("difFinitOut")),

            tabItem("Derivacion2",
                    h1("Diferencias Finitas Gradiente"),
                    box(textInput("difFinEcu_2", "Ingrese la Ecuación"),
                        textInput("valorX_2", "x"),
                        textInput("valorY_2", "y"),
                        textInput("valorH_2", "h")),
                    actionButton("diferFinEval_2", "Evaluar Gradiente"),
                    tableOutput("difFinitOut_2")),
            tabItem("newton_method",
                    h1("Metodo de Newton-Raphson"),
                    box(textInput("ec_nwt", "Ingrese la Ecuación"),
                        textInput("x0_nwt", "X0"),
                        textInput("k_max_nwt", "Kmax"),
                        textInput("e_nwt", "Error")),
                    actionButton("nwtRaphsonSolver", "Newton-Raphson Solver"),
                    tableOutput("output_nwt_df")),
            tabItem("bisection_method",
                    h1("Metodo de Biseccion"),
                    box(textInput("ec_bisct", "Ingrese la Ecuación"),
                        textInput("a_bisct", "a"),
                        textInput("b_bisct", "b"),
                        textInput("k_max_bisct", "Kmax"),
                        textInput("e_bisct", "Error")),
                    actionButton("bisctSolver", "Bisection Solver"),
                    tableOutput("output_bisct_df")),
            tabItem("GD_Q",
                    h1("Metodo de GD (Q)"),
                    box(textInput("q_gd", "Ingrese Q en el formato [[a1,a2],[a3,a4]]"),
                        textInput("c_gd", "Ingrese C en el formato [c1,c2]"),
                        textInput("x0_gd", "Ingrese X0 en el formato [a,b]"),
                        textInput("e_gd", "Ingrese Tolerancia"),
                        textInput("n_gd", "N (maximo iteraciones)"),
                        textInput("lr_gd","LR (Solo Metodo Constante)"),
                        textInput("lr_method_gd", "Metodo de Learning Rate")),
                    actionButton("gradieDescentQSolver", "Quadratic Function GD Solver"),
                    tableOutput("output_gd_q_df")),
            tabItem("GD_Rosenbrock",
                    h1("Metodo de GD (Rosenbrock)"),
                    box(textInput("xo_gd2", "Ingrese Q en el formato [a1,a2]"),
                        textInput("e_gd2", "Ingrese Tolerancia"),
                        textInput("n_gd2", "N (maximo iteraciones)"),
                        textInput("lr_gd2", "Ingrese LR")),
                    actionButton("gradieDescentRosenbrockSolver", "Rosenbrock Function GD Solver"),
                    tableOutput("output_gd_rosenbrock_df")),
            tabItem("gen_data",
                    h1("Generate Data"),
                    box(textInput("d", "Ingrese numero columnas (d)"),
                        textInput("n", "Ingrese numero observaciones (n)"),
                        textInput("path", "Ingrese ruta")),
                    actionButton("gen_data_btn", "Generate"),
                    textOutput("output_gen_data")),
            tabItem("Closed_Solution",
                    h1("Solucion Cerrada"),
                    box(
                        textInput("path_closed", "Ingrese ruta")),
                    actionButton("closedSolutionSolver", "Closed Solution"),
                    tableOutput("output_closed_solution")),
            tabItem("GD_Rosenbrock2",
                    h1("Metodo de GD con Backtracking (Rosenbrock)"),
                    box(textInput("xo_rb_bt_gd", "Ingrese X0 en el formato [a1,a2]"),
                        textInput("e_rb_bt_gd", "Ingrese Tolerancia"),
                        textInput("n_rb_bt_gd", "N (maximo iteraciones)"),
                        textInput("lr_rb_bt_gd", "Ingrese LR")),
                    actionButton("rosenbrack_backtracking_gdSolver", "Rosenbrock Function GD Solver"),
                    tableOutput("out_rosenbrock_backtracking_gd")),
            tabItem("Newton_Rosenbrock",
                    h1("Metodo de Newton con Backtracking (Rosenbrock)"),
                    box(textInput("xo_rb_bt_nwtn", "Ingrese X0 en el formato [a1,a2]"),
                        textInput("e_rb_bt_nwtn", "Ingrese Tolerancia"),
                        textInput("n_rb_bt_nwtn", "N (maximo iteraciones)"),
                        textInput("lr_rb_bt_nwtn", "Ingrese LR")),
                    actionButton("rosenbrack_backtracking_newtonSolver", "Rosenbrock Function Newton Solver"),
                    tableOutput("out_rosenbrock_backtracking_nwtn")),
            tabItem("bgd_tab",
                    h1("Metodo de GD"),
                    box(
                        textInput("e_bgd", "Ingrese Tolerancia"),
                        textInput("n_bgd", "N (maximo iteraciones)"),
                        textInput("lr_bgd", "Ingrese LR"),
                        textInput("path_bgd", "Ingrese Ruta")),
                    actionButton("batchGDSolver", "Batch GD Solver"),
                    tableOutput("out_batch_gd")),
            tabItem("sgd_tab",
                    h1("Metodo SGD"),
                    box(
                        textInput("e_sgd", "Ingrese Tolerancia"),
                        textInput("n_sgd", "N (maximo iteraciones)"),
                        textInput("lr_sgd", "Ingrese LR"),
                        textInput("path_sgd", "Ingrese Ruta")),
                    actionButton("sgdSolver", "SGD Solver"),
                    tableOutput("out_sgd")),
            tabItem("mbgd_tab",
                    h1("Metodo Mini Batch GD"),
                    box(
                        textInput("e_mbgd", "Ingrese Tolerancia"),
                        textInput("n_mbgd", "N (maximo iteraciones)"),
                        textInput("lr_mbgd", "Ingrese LR"),
                        textInput("bs_mbgd", "Ingrese Batch Size"),
                        textInput("path_mbgd", "Ingrese Ruta")),
                    actionButton("mbgdSolver", "MBGD Solver"),
                    tableOutput("out_mbgd"))
        )
    )
)

