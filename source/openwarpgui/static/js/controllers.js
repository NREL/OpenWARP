'use strict';
/*jshint -W097*/
/*jshint strict:false*/
/*global angular, document, setTimeout*/

// controllers
angular.module('openwarp.controllers', [])
//init controller
	.controller('init', ['$scope', '$rootScope', '$window', function ($scope, $rootScope, $window) {

	} ])
// meshing controller
    .controller('meshing', ['$scope', '$rootScope', '$window', 'localStorageService', function ($scope, $rootScope, $window, localStorageService) {

        //retrieve text input data from local storage
        $("input.text").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);

            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
        });

        //retrieve radio/checkbox data from local storage
        $("input.radio, input.checkbox").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
            if ($scope[id] === 'true') {
                $(this).prop("checked", true);
            } else {
                $(this).prop("checked", false);
            }
        });

        //save input value into local storage after change
        $("input.text").change(function () {
            if (!$(this).hasClass("ban-storage")) {
                localStorageService.set($(this).attr("id"), $(this).val());
            }
        });

        //save radio/checkbox value
        $("input.radio, input.checkbox").change(function () {
            var wrapper = $(this).parents(".radio-wrapper").eq(0);
            $("input.radio, input.checkbox", wrapper).each(function () {
                var radioStaus;
                if ($(this).prop("checked")) {
                    radioStaus = "true";
                } else {
                    radioStaus = "false";
                }
                localStorageService.set($(this).attr("id"), radioStaus);
            });
        });

        // Invoked when the MESHING EXECUTE button was clicked.
        $("#meshing_execute").click(function () {
            // Loading label to prevent editing or clicking operations
            openWARP.app.show_loading();
            // Collect parameters
            var general_keys = ["infile", "outfile", "maxh", "minh", "fineness",
                            "grading", "tolerance"];
            var data_to_post = {};
            var invalidate = false;
            $.each(general_keys, function (_itr, item) {
                // Check if the value is accepted locally by input
                if ($("#meshing_" + item).hasClass('input-error')) {
                    invalidate = true;
                    return false;
                }
                data_to_post[item] = $("#meshing_" + item).val();
            });
            if (invalidate) {
                openWARP.app.message_box("Input error:", "Check the red input boxes to see the invalid values.");
                return;
            }
            data_to_post["usetolerance"] = $('input[name="Usetolerance"]:checked').val();
            if (!(data_to_post.usetolerance === "1")) {
                // Assign any value to "tolerance"
                data_to_post["tolerance"] = "0";
            }
            // Validate parameters
            if (!data_to_post.infile || data_to_post.infile === "") {
                openWARP.app.message_box("Input error:", "Input file must not be empty.");
                return;
            }
            if (!data_to_post.outfile || data_to_post.outfile === "") {
                openWARP.app.message_box("Input error:", "Output file must not be empty.");
                return;
            }
            if (!data_to_post.maxh || data_to_post.maxh === "") {
                openWARP.app.message_box("Input error:", "Maxh must not be empty.");
                return;
            }
            if (!data_to_post.minh || data_to_post.minh === "") {
                openWARP.app.message_box("Input error:", "Minh must not be empty.");
                return;
            }
            if (!data_to_post.fineness || data_to_post.fineness === "") {
                openWARP.app.message_box("Input error:", "Finess file must not be empty.");
                return;
            }
            if (!data_to_post.grading || data_to_post.grading === "") {
                openWARP.app.message_box("Input error:", "Grading must not be empty.");
                return;
            }
            if (!data_to_post.usetolerance || data_to_post.usetolerance === "") {
                openWARP.app.message_box("Input error:", "Usetolerance radios must be selected.");
                return;
            }
            if (data_to_post.usetolerance === "1" && (!data_to_post.tolerance || data_to_post.usetolerance === "")) {
                openWARP.app.message_box("Input error:", "Tolerance must not be empty.");
                return;
            }
            if (parseFloat(data_to_post.maxh) < parseFloat(data_to_post.minh)) {
                openWARP.app.message_box("Input error:", "Minh must not be larger than Maxh.");
                return;
            }

            // AJAX send request
            $.ajax({
                url: "/generate_mesh",
                type: "POST",
                dataType: "json",
                data: data_to_post,
                success: function (data) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // extract log and display in overlay dialog
                    var logContent = data.log;
                    openWARP.app.message_box("Running log in server side:", logContent);
                },
                error: function (request, status, error) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        });

        // Quit application
        $("#application_quit").click(function () {
            openWARP.app.quit_app();
        });

        //reset meshing page 
        openWARP.app.screenReset("meshing");

    } ])
// simulation controller
    .controller('simulation', ['$scope', '$rootScope', '$window', 'localStorageService', function ($scope, $rootScope, $window, localStorageService) {

        //retrieve text input data from local storage
        $("input.text").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
        });

        //retrieve radio/checkbox data from local storage
        $("input.radio, input.checkbox").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
            if ($scope[id] === 'true') {
                $(this).prop("checked", true);
            } else {
                $(this).prop("checked", false);
            }
        });

        //save input value into local storage after change
        $("input.text").change(function () {
            if (!$(this).hasClass("ban-storage")) {
                localStorageService.set($(this).attr("id"), $(this).val());
            }
        });

        //save radio/checkbox value
        $("input.radio, input.checkbox").change(function () {
            var wrapper = $(this).parents(".radio-wrapper").eq(0);
            $("input.radio, input.checkbox", wrapper).each(function () {
                var radioStaus;
                if ($(this).prop("checked")) {
                    radioStaus = "true";
                } else {
                    radioStaus = "false";
                }
                localStorageService.set($(this).attr("id"), radioStaus);
            });
        });

        //Set default to 5 FLOATIING BODY 
        if ($scope.simulation_body_num === null) {
            $("#simulation_body_num").val(5).trigger("change");
            $scope.simulation_body_num = 5;
        }

        //show saved body number tabs
        var bodyNum = parseInt($scope.simulation_body_num, 10);
        $(".tab-bodies ul li").each(function (index, element) {
            if (index < bodyNum) {
                $(this).removeClass("hide");
            } else {
                $(this).addClass("hide");
            }
        });
        if (bodyNum === 10) {
            $(".js-add-body").hide();
        }

        //floating bodies tab click 
        $(".tab-bodies ul li a").click(function () {
            if (!$(this).parent().hasClass("active")) {
                $(".tab-bodies ul li").removeClass("active");
                $(this).parent().addClass("active");
                $(".tab-body").addClass("hide");
                var idx = $(this).parent().index();
                var tTab = $(".tab-body").eq(idx);
                tTab.removeClass("hide");
                $("input.text-num", tTab).trigger("keyup");
                $("input.text-filename", tTab).trigger("blur");
                $(".input-wrapper .unit", tTab).each(function () {
                    if ($(this).prev("input.text").val() === "") {
                        $(this).show();
                    } else {
                        $(this).hide();
                    }
                });
            }
        });

        //add floating body
        $(".js-add-body").click(function () {
            bodyNum = parseInt($scope.simulation_body_num, 10);
            $scope.simulation_body_num = bodyNum + 1;
            $(".tab-bodies ul li").eq(bodyNum).removeClass("hide").find("a").trigger("click");
            $("#simulation_body_num").val(bodyNum + 1).trigger("change");
            if (bodyNum === 9) {
                $(this).hide();
            }
        });

        // Invoked when SIMULATION EXECUTE button is clicked.
        $("#simulation_execute").click(function () {
            openWARP.app.show_loading();
            // Checking parameters locally
            var invalidate = false;
            $('input').each(function () {
                if ($(this).hasClass('input-error')) {
                    openWARP.app.message_box("Input error:", "The field " + $(this).attr('id') + " is invalid.");
                    invalidate = true;
                    return false;
                }
            });
            if (invalidate) {
                // Loading label is already hidden by the "Input error" warning.
                return;
            }

            // Construct general parameters
            var general_keys = ["rho", "g", "depth", "xeff", "yeff",
                            "wave_frequencies", "min_wave_frequencies", "max_wave_frequencies",
                            "wave_directions", "max_wave_direction", "min_wave_directions",
                            "indiq_solver", "ires", "tol_gmres", "max_iterations", "save_potential",
                            "green_tabulation_numx", "green_tabulation_numz", "green_tabulation_simpson_npoints",
                            "use_ode_influence_coefficients", "use_higher_order", "num_panel_higher_order",
                            "b_spline_order", "use_dipoles_implementation", "thin_panels",
                            "compute_drift_forces", "compute_yaw_moment", "remove_irregular_frequencies"];
            var data_to_post = {};
            $.each(general_keys, function (_itr, item) {
                data_to_post[item] = $("#simulation_" + item).val();
                if (!data_to_post[item] || data_to_post[item] === "") {
                    openWARP.app.message_box("Input error:", item + " must not be empty.");
                    invalidate = true;
                    return false;
                }
            });
            if (invalidate) {
                // Loading label is already hidden by the "Input error" warning.
                return;
            }
            // Check minimum and maximum
            if (parseFloat(data_to_post.min_wave_frequencies) > parseFloat(data_to_post.max_wave_frequencies)) {
                openWARP.app.message_box("Input error:", "min_wave_frequencies must not be larger than max_wave_frequencies.");
                return;
            }
            if (parseFloat(data_to_post.min_wave_directions) > parseFloat(data_to_post.max_wave_direction)) {
                openWARP.app.message_box("Input error:", "min_wave_directions must not be larger than max_wave_directions.");
                return;
            }

            // Construct bodies parameters
            var bodies = [];
            var body_general_keys = ["mesh_file", "points", "panels", "additional_info_lines"];
            for (var i = 1; i <= $('.tab-body').length; i++) {
                var body = {};
                // Check if post this floating body or not
                if ($("#simulation_b" + i + "_mesh_file").val() === "") {
                    continue;
                }
                // general body keys
                $.each(body_general_keys, function (_itr, item) {
                    body[item] = $("#simulation_b" + i + "_" + item).val();
                    if (!body[item] || body[item] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": " + item + " must not be empty.");
                        invalidate = true;
                        return false;
                    }
                });
                if (invalidate) {
                    return;
                }

                // DOFs and forces
                var num_dof = 0;
                var num_forces = 0;
                var prefix_line;
                // surge
                prefix_line = "#simulation_b" + i + "_surge";
                if ($(prefix_line).is(':checked')) {
                    body["surge"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 0);
                    if (body["surge"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": surge vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["surge"] = null;
                }
                // sway
                prefix_line = "#simulation_b" + i + "_sway";
                if ($(prefix_line).is(':checked')) {
                    body["sway"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 1);
                    if (body["sway"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": sway vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["sway"] = null;
                }
                // heave
                prefix_line = "#simulation_b" + i + "_heave";
                if ($(prefix_line).is(':checked')) {
                    body["heave"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 2);
                    if (body["heave"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": heave vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["heave"] = null;
                }
                // roll_about_cdg
                prefix_line = "#simulation_b" + i + "_roll_about_cdg";
                if ($(prefix_line).is(':checked')) {
                    body["roll_about_cdg"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 0);
                    if (body["roll_about_cdg"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": roll_about_cdg vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["roll_about_cdg"] = null;
                }
                // pitch_about_cdg
                prefix_line = "#simulation_b" + i + "_pitch_about_cdg";
                if ($(prefix_line).is(':checked')) {
                    body["pitch_about_cdg"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 1);
                    if (body["pitch_about_cdg"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": pitch_about_cdg vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["pitch_about_cdg"] = null;
                }
                // yaw_about_cdg
                prefix_line = "#simulation_b" + i + "_yaw_about_cdg";
                if ($(prefix_line).is(':checked')) {
                    body["yaw_about_cdg"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 2);
                    if (body["yaw_about_cdg"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": yaw_about_cdg vector must be all set properly.");
                        return;
                    }
                    num_dof++;
                } else {
                    body["yaw_about_cdg"] = null;
                }
                // force_in_x_direction
                prefix_line = "#simulation_b" + i + "_force_in_x_direction";
                if ($(prefix_line).is(':checked')) {
                    body["force_in_x_direction"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 0);
                    if (body["force_in_x_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": force_in_x_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["force_in_x_direction"] = null;
                }
                // force_in_y_direction
                prefix_line = "#simulation_b" + i + "_force_in_y_direction";
                if ($(prefix_line).is(':checked')) {
                    body["force_in_y_direction"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 1);
                    if (body["force_in_y_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": force_in_y_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["force_in_y_direction"] = null;
                }
                // force_in_z_direction
                prefix_line = "#simulation_b" + i + "_force_in_z_direction";
                if ($(prefix_line).is(':checked')) {
                    body["force_in_z_direction"] = openWARP.app.construct3dVectorParameter(true, prefix_line, 2);
                    if (body["force_in_z_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": force_in_z_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["force_in_z_direction"] = null;
                }
                // moment_cdg_force_in_x_direction
                prefix_line = "#simulation_b" + i + "_moment_cdg_force_in_x_direction";
                if ($(prefix_line).is(':checked')) {
                    body["moment_cdg_force_in_x_direction"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 0);
                    if (body["moment_cdg_force_in_x_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": moment_cdg_force_in_x_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["moment_cdg_force_in_x_direction"] = null;
                }
                // moment_cdg_force_in_y_direction
                prefix_line = "#simulation_b" + i + "_moment_cdg_force_in_y_direction";
                if ($(prefix_line).is(':checked')) {
                    body["moment_cdg_force_in_y_direction"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 1);
                    if (body["moment_cdg_force_in_y_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": moment_cdg_force_in_y_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["moment_cdg_force_in_y_direction"] = null;
                }
                // moment_cdg_force_in_z_direction
                prefix_line = "#simulation_b" + i + "_moment_cdg_force_in_z_direction";
                if ($(prefix_line).is(':checked')) {
                    body["moment_cdg_force_in_z_direction"] = openWARP.app.construct3dVectorParameter(false, prefix_line, 2);
                    if (body["moment_cdg_force_in_z_direction"] === "") {
                        openWARP.app.message_box("Input error:", "Body" + i + ": moment_cdg_force_in_z_direction vector must be all set properly.");
                        return;
                    }
                    num_forces++;
                } else {
                    body["moment_cdg_force_in_z_direction"] = null;
                }

                body["degrees_of_freedom"] = "" + num_dof;
                body["resulting_generalised_forces"] = "" + num_forces;

                // Append this floating body
                bodies.push(body);
            }
            // Merge parameters
            data_to_post["floating_bodies"] = bodies;

            // AJAX send request
            $.ajax({
                url: "/simulate",
                type: "POST",
                dataType: "json",
                data: {
                    json_str: JSON.stringify(data_to_post)
                },
                success: function (data) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // extract log and display in overlay dialog
                    var logContent = data.log;
                    openWARP.app.message_box("Running log in server side:", logContent);
                },
                error: function (request, status, error) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        });

        // Quit application
        $("#application_quit").click(function () {
            openWARP.app.quit_app();
        });

        //reset simulation page 
        openWARP.app.screenReset("simulation");
    } ])
// postprocessing controller
    .controller('postprocessing', ['$scope', '$rootScope', '$window', 'localStorageService', function ($scope, $rootScope, $window, localStorageService) {

        //retrieve text input data from local storage
        $("input.text").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
        });

        //retrieve radio/checkbox data from local storage
        $("input.radio, input.checkbox").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
            if ($scope[id] === 'true') {
                $(this).prop("checked", true);
            } else {
                $(this).prop("checked", false);
            }
        });

        //save input value into local storage after change
        $("input.text").change(function () {
            if (!$(this).hasClass("ban-storage")) {
                localStorageService.set($(this).attr("id"), $(this).val());
            }
        });

        //save radio/checkbox value
        $("input.radio, input.checkbox").change(function () {
            var wrapper = $(this).parents(".radio-wrapper").eq(0);
            $("input.radio, input.checkbox", wrapper).each(function () {
                var radioStaus;
                if ($(this).prop("checked")) {
                    radioStaus = "true";
                } else {
                    radioStaus = "false";
                }
                localStorageService.set($(this).attr("id"), radioStaus);
            });
        });

        // Save as Tecplot button clicked
        $("#postprocess_execute").click(function () {
            openWARP.app.show_loading();
            // Checking parameters locally
            var invalidate = false;
            $('input').each(function () {
                if ($(this).hasClass('input-error')) {
                    openWARP.app.message_box("Input error:", "The field " + $(this).attr('id') + " is invalid.");
                    invalidate = true;
                    return false;
                }
            });
            if (invalidate) {
                // Loading label is already hidden by the "Input error" warning.
                return;
            }

            // Construct parameters
            var data_to_post = {};
            // irf
            var use_irf = $('input[name="irf"]:checked').val();
            if (!use_irf || use_irf === "") {
                openWARP.app.message_box("Input error:", "Use irf radios must be selected.");
                return;
            }
            if (use_irf === "1") {
                var irf = [];
                // Use irf
                irf[0] = use_irf;
                // time step
                var time_step = $("#postprocessing_time_step").val();
                if (!time_step || time_step === "") {
                    openWARP.app.message_box("Input error:", "time step must not be empty.");
                    return;
                }
                if (parseFloat(time_step) < 0) {
                    openWARP.app.message_box("Input error:", "time step must not be negative.");
                    return;
                }
                irf[1] = time_step;
                // duration
                var duration = $("#postprocessing_duration").val();
                if (!duration || duration === "") {
                    openWARP.app.message_box("Input error:", "duration must not be empty.");
                    return;
                }
                if (parseFloat(duration) < 0) {
                    openWARP.app.message_box("Input error:", "duration must not be negative.");
                    return;
                }
                irf[2] = duration;
                data_to_post["irf"] = irf;
            } else {
                data_to_post["irf"] = [use_irf];
            }
            // show_pressure
            data_to_post["show_pressure"] = $('input[name="showPressure"]:checked').val();
            if (!data_to_post["show_pressure"] || data_to_post["show_pressure"] === "") {
                openWARP.app.message_box("Input error:", "Show_pressure radios must be selected.");
                return;
            }
            // Some general keys
            var general_keys = ["kochin_function", "min_angle", "max_angle", "number_points_x",
                                "number_points_y", "dimensions_x", "dimensions_y"];
            var general_collector = {};
            $.each(general_keys, function (_itr, item) {
                general_collector[item] = $("#postprocessing_" + item).val();
                if (!general_collector[item] || general_collector[item] === "") {
                    openWARP.app.message_box("Input error:", item + " must not be empty.");
                    invalidate = true;
                    return false;
                }
            });
            if (invalidate) {
                // Loading label is already hidden by the "Input error" warning.
                return;
            }
            // Check minimum and maximum
            if (parseFloat(general_collector.min_angle) > parseFloat(general_collector.max_angle)) {
                openWARP.app.message_box("Input error:", "min_angle must not be larger than max_angle.");
                return;
            }
            // Assign arrays
            if (general_collector.kochin_function === "0") {
                data_to_post["kochin_function"] = [general_collector.kochin_function];
            } else {
                data_to_post["kochin_function"] = [general_collector.kochin_function,
                                                   general_collector.min_angle,
                                                   general_collector.max_angle];
            }
            data_to_post["free_surface_elevation"] = [general_collector.number_points_x,
                                                        general_collector.number_points_y,
                                                        general_collector.dimensions_x,
                                                        general_collector.dimensions_y];

            // Check positive and integer
            if (!openWARP.app.isNonNegativeInteger(general_collector.kochin_function)) {
                openWARP.app.message_box("Input error:", "kochin_function must be non-negative integer.");
                return;
            }
            if (!openWARP.app.isNonNegativeInteger(general_collector.number_points_x)) {
                openWARP.app.message_box("Input error:", "number_points_x must be non-negative integer.");
                return;
            }
            if (!openWARP.app.isNonNegativeInteger(general_collector.number_points_y)) {
                openWARP.app.message_box("Input error:", "number_points_y must be non-negative integer.");
                return;
            }
            if (!openWARP.app.isNonNegativeInteger(general_collector.dimensions_x)) {
                openWARP.app.message_box("Input error:", "dimensions_x must be non-negative integer.");
                return;
            }
            if (!openWARP.app.isNonNegativeInteger(general_collector.dimensions_y)) {
                openWARP.app.message_box("Input error:", "dimensions_y must be non-negative integer.");
                return;
            }
            if (!openWARP.app.isNonNegativeInteger(general_collector.kochin_function)) {
                openWARP.app.message_box("Input error:", "kochin_function must be non-negative integer.");
                return;
            }

            // Send AJAX request
            $.ajax({
                url: "/postprocess",
                type: "POST",
                dataType: "json",
                data: {
                    json_str: JSON.stringify(data_to_post)
                },
                success: function (data) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // extract log and display in overlay dialog
                    var logContent = data.log;
                    openWARP.app.message_box("Running log in server side:", logContent);
                },
                error: function (request, status, error) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        });

        // Visualize the simulation result
        $("#visualize_simulation").click(function () {
            openWARP.app.show_loading();
            // Send AJAX request
            $.ajax({
                url: "/visualize",
                type: "POST",
                success: function (data) {
                    openWARP.app.message_box("Visulization:", "Launched ...");
                },
                error: function (request, status, error) {
                    // Notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        });

        // Quit application
        $("#application_quit").click(function () {
            openWARP.app.quit_app();
        });

        //reset postprocessing page 
        openWARP.app.screenReset("postprocessing");
    } ])
	// configuration controller
    .controller('configuration', ['$scope', '$rootScope', '$window', 'localStorageService', function ($scope, $rootScope, $window, localStorageService) {

        //retrieve text input data from local storage
        $("input.text").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);

            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
        });

        //retrieve radio/checkbox data from local storage
        $("input.radio, input.checkbox").each(function () {
            var id = $(this).attr("id");
            $scope[id] = localStorageService.get(id);
            if ($(this)[0].hasAttribute("data-default-value")) {
                if (($scope[id] === null) || ($scope[id] === undefined) || ($scope[id] === "")) {
                    $scope[id] = $(this).attr("data-default-value"); ;
                }
            }
            if ($scope[id] === 'true') {
                $(this).prop("checked", true);
            } else {
                $(this).prop("checked", false);
            }
        });

        //save input value into local storage after change
        $("input.text").change(function () {
            if (!$(this).hasClass("ban-storage")) {
                localStorageService.set($(this).attr("id"), $(this).val());
            }
        });

        //save radio/checkbox value
        $("input.radio, input.checkbox").change(function () {
            var wrapper = $(this).parents(".radio-wrapper").eq(0);
            $("input.radio, input.checkbox", wrapper).each(function () {
                var radioStaus;
                if ($(this).prop("checked")) {
                    radioStaus = "true";
                } else {
                    radioStaus = "false";
                }
                localStorageService.set($(this).attr("id"), radioStaus);
            });
        });

        // Invoked when the MESHING EXECUTE button was clicked.
        $("#apply_configuration_changes").click(function () {
            // Loading label to prevent editing or clicking operations
            openWARP.app.show_loading();
            // Collect parameters

            var data_to_post = {}
            // Logging level
            data_to_post["logging_level"] = $('input[name="logging_level"]:checked').val();
            if (!(data_to_post.logging_level === "20")) {
                data_to_post["logging_level"] = "10";
            }

	        // Clear log flag
	        data_to_post["clear_log_flag"] = $('input[name="clear_log_flag"]:checked').val();
            if (!(data_to_post["clear_log_flag"] === "true")) {
                data_to_post["clear_log_flag"] = "false";
            }

            // AJAX send request
            $.ajax({
                url: "/apply_configuration",
                type: "POST",
                dataType: "json",
                data: data_to_post,
                success: function (data) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // extract log and display in overlay dialog
                    var logContent = data.log;
                    openWARP.app.message_box("Running log in server side:", logContent);
                },
                error: function (request, status, error) {
                    // Explicitly hide loading dialog
                    openWARP.app.hide_loading();
                    // notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        });

        // Quit application
        $("#application_quit").click(function () {
            openWARP.app.quit_app();
        });

        //reset meshing page
        openWARP.app.screenReset("configuration");

    } ])
	