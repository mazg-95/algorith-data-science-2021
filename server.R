
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

    # ##################################################

    generateData <- eventReactive(input$gen_data_btn, {
        d_str <- input$d[1]
        n_str <- input$n[1]
        path_str <- input$path[1]
        outs <- generate_and_save_data(d_str, n_str, path_str)
        "Success"
    })
    

    closed_solutionCalc <- eventReactive(input$closedSolutionSolver, {
        path_str <- input$path_closed[1]
        outs <- closed_solution_solver(path_str)
        outs
    })

    batchGDCalc <- eventReactive(input$batchGDSolver, {
        e_str <- input$e_bgd[1]
        n_str <- input$n_bgd[1]
        lr_str <- input$lr_bgd[1]
        path <- input$path_bgd[1]
        gd_solver(e_str, n_str, lr_str, path=path)}
    )

    sgdCalc <- eventReactive(input$sgdSolver, {
        e_str <- input$e_sgd[1]
        n_str <- input$n_sgd[1]
        lr_str <- input$lr_sgd[1]
        path <- input$path_sgd[1]
        browser()
        gd_solver(e_str, n_str, lr_str, path=path, method="sgd")}
    )

    mbgdCalc <- eventReactive(input$mbgdSolver, {
        e_str <- input$e_mbgd[1]
        n_str <- input$n_mbgd[1]
        lr_str <- input$lr_mbgd[1]
        batch_size_str <- input$bs_mbgd[1]
        path <- input$path_mbgd[1]
        gd_solver(e_str, n_str, lr_str, batch_size_str=batch_size_str, path=path, method="mbgd")}
    )

    rosenbrock_backtracking_gdCalc  <- eventReactive(input$rosenbrack_backtracking_gdSolver, {
        x0_str <- input$xo_rb_bt_gd[1]
        e_str <- input$e_rb_bt_gd[1]
        n_str <- input$n_rb_bt_gd[1]
        lr_str <- input$lr_rb_bt_gd[1]
        ronsenbrock_gd_backtracking_solver(x0_str, e_str, n_str, lr_str)}
    )

    rosenbrock_backtracking_nwtnCalc  <- eventReactive(input$rosenbrack_backtracking_newtonSolver, {
        x0_str <- input$xo_rb_bt_nwtn[1]
        e_str <- input$e_rb_bt_nwtn[1]
        n_str <- input$n_rb_bt_nwtn[1]
        lr_str <- input$lr_rb_bt_nwtn[1]
        ronsenbrock_newton_backtracking_solver(x0_str, e_str, n_str, lr_str)}
    )


    output$output_gen_data <- renderText({
        generateData()
    })
    
    output$output_closed_solution <- renderTable({
        closed_solutionCalc()
    })

    output$out_batch_gd <- renderTable({
        batchGDCalc()
    })

    output$out_sgd <- renderTable({
        sgdCalc()
    })

    output$out_mbgd <- renderTable({
        mbgdCalc()
    })

    output$out_rosenbrock_backtracking_gd <- renderTable({
        rosenbrock_backtracking_gdCalc()
    })

    output$out_rosenbrock_backtracking_nwtn <- renderTable({
        rosenbrock_backtracking_nwtnCalc()
    })

    # ##########################################

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
