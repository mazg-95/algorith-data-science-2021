
library(shiny)
library(reticulate)

source_python("algoritmos.py")

#tableOut, soluc = newtonSolverX(-5, "2x^5 - 3", 0.0001)


shinyServer(function(input, output) {
    
    #Evento y evaluación de metodo de newton para ceros
    newtonCalculate<-eventReactive(input$nwtSolver, {
        inputEcStr<-input$ecuacion[1]
        initVal<-input$initVal[1]
        error<-input$Error[1]
        outs<-newtonSolverX(initVal, inputEcStr, error)
        outs
    })
    
    #Evento y evaluación de diferencias finitas
    diferFinitCalculate<-eventReactive(input$diferFinEval, {
        inputEcStr<-input$difFinEcu[1]
        valX<-input$valorX[1]
        h<-input$valorH[1]
        outs <- evaluate_finite_difference(inputEcStr, valX, h)
        outs
    })
    
    diferFinitCalculate2<-eventReactive(input$diferFinEval_2, {
        inputEcStr<-input$difFinEcu_2[1]
        valX<-input$valorX_2[1]
        valY <- input$valorY_2[1]
        h<-input$valorH_2[1]
        outs<-evaluate_finite_difference_2(inputEcStr, valX, valY, h)
        outs
    })
    
    newtonRaphsonCalc<-eventReactive(input$nwtRaphsonSolver, {
        inputEcStr<-input$ec_nwt[1]
        x0<-input$x0_nwt[1]
        k_max <- input$k_max_nwt[1]
        error<-input$e_nwt[1]
        outs<-newton_raphson_solver(inputEcStr, x0, k_max, error)
        outs
    })
    
    bisectionCalc<-eventReactive(input$bisctSolver, {
        inputEcStr<-input$ec_bisct[1]
        a<-input$a_bisct[1]
        b<-input$b_bisct[1]
        k_max <- input$k_max_bisct[1]
        error<-input$e_bisct[1]
        outs<-bisection_solver(inputEcStr, a, b, k_max, error)
        outs
    })

    q_gdCalc<-eventReactive(input$gradieDescentQSolver, {
        q_str<-input$q_gd[1]
        c_str<-input$c_gd[1]
        e_str<-input$e_gd[1]
        n_str<-input$n_gd[1]
        lr_str<-input$lr_gd[1]
        lr_method<-input$lr_method_gd[1]
        x0_str<-input$x0_gd[1]
        outs<-quadratic_gd_solver(q_str, c_str, x0_str, e_str, n_str, lr_str, lr_method)
        outs
    })

    rosenbrock_gdCalc<-eventReactive(input$gradieDescentRosenbrockSolver, {
        x0_str<-input$xo_gd2[1]
        e_str<-input$e_gd2[1]
        n_str<-input$n_gd2[1]
        lr_str<-input$lr_gd2[1]
        outs<-ronsenbrock_gd_solver(x0_str, e_str, n_str, lr_str)
        outs
    })

    
    #REnder metodo de Newton
    output$salidaTabla<-renderTable({
        newtonCalculate()
    })
    
    #Render Diferncias Finitas 1
    output$difFinitOut<-renderTable({
        diferFinitCalculate()
    }, digits = 4)
    
    #Render Diferncias Finitas 2
    output$difFinitOut_2<-renderTable({
        diferFinitCalculate2()
    }, digitls = 4)

    output$output_nwt_df<-renderTable({
        newtonRaphsonCalc()
    }, digits = 5)
    output$output_bisct_df<-renderTable({
        bisectionCalc()
    }, digits = 5)

    output$output_gd_q_df<-renderTable({
        q_gdCalc()
    })

    output$output_gd_rosenbrock_df<-renderTable({
        rosenbrock_gdCalc()
    })
    
})
