program main

    use test_metric_mod,       only : test_metric, test_metric_vert
    use test_metric_class_mod, only : test_metric_class
    use parcomm_mod,           only : init_global_parallel_enviroment, &
                                      deinit_global_parallel_enviroment

    call init_global_parallel_enviroment()

    call test_metric_class("cube","ecs")
    call test_metric_class("cube","shallow_atmosphere_metric")
    call test_metric()
    call test_metric_vert()

    call deinit_global_parallel_enviroment()

end
