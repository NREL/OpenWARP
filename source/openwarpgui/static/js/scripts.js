// Use this variable to set up the common and page specific functions.
var openWARP = {
    // All pages
    common: {
        init: function () {
            // JavaScript to be fired on all pages
        }
    },
    landing: {
        init: function () {
            // JavaScript to be fired on landing pages
            // Quit application
            $("#application_quit").click(function () {
                openWARP.app.quit_app();
            });
        }
    },
    app: {
        init: function () {
            // JavaScript to be fired on app pages
            //set page min-height according to left side navigation height
            $("html, body").css({
                "min-height": parseInt($(".app-nav").height(), 10) + 183
            });

            //set height after resize
            $(window).resize(function () {
                $(".openwarp-app, .app-screen").height("");
                var h = Math.max((parseInt($(".app-screen").height(), 10)), (parseInt($(".openwarp-app").height(), 10) - 57));
                $(".openwarp-app").height(h + 57);
                $(".app-screen").height(h);
            });

            // Setup the jQuery dialogs
            $("#single_dialog").dialog({
                autoOpen: false,
                title: 'OpenWARP information',
                modal: true,
                resizable: false,
                width: 600,
                closeOnEscape: false
            });
        },
        heightReset: function () {
            $('html, body').animate({ scrollTop: 0 });
            $(".openwarp-app, .app-screen").height("");
            var h = Math.max((parseInt($(".app-screen").height(), 10)), (parseInt($(".openwarp-app").height(), 10) - 57));
            $(".openwarp-app").height(h + 57);
            $(".app-screen").height(h);
        },
        screenReset: function (navItem) {
            $(".app-nav li").removeClass("on");
            $(".app-nav li.link-" + navItem).addClass("on");
            openWARP.app.heightReset();

            //init radio
            $(".radio-wrapper input.radio, .radio-wrapper input.checkbox").each(function () {
                var forId = $(this).attr('id');
                if ($(this).prop("checked")) {

                    $(".input-dot[for='" + forId + "']").addClass("input-dot-checked");
                }
            });

            //disable input according to radio check
            $(".js-disable-fieldset").each(function () {
                if ($(this).prop("checked")) {
                    var tId = $(this).attr('data-ctrl-target');
                    var wrapper = $("." + tId);
                    $("input.text", wrapper).val("").prop("disabled", true).removeClass("input-error"); ;
                    $("fieldset", wrapper).addClass("fieldset-disabled");
                }
            });
            //disable input according to radio check
            $(".js-toggle-fieldset").each(function () {
                if (!$(this).prop("checked")) {
                    var tId = $(this).attr('data-ctrl-target');
                    var wrapper = $("." + tId);
                    $("input.text", wrapper).val("").prop("disabled", true).removeClass("input-error"); ;
                    wrapper.addClass("fieldset-disabled");
                }
            });

            //prevent help dropdown overlap
            $(".screen-sub-title").each(function (index, element) {
                $(this).css({
                    "z-index": (30 - index)
                });
            });

            //Show/hide Help Dropdown
            $(".help-wrapper .help-ico, .help-wrapper .close-ico").click(function () {
                var wrapper = $(this).parents(".help-wrapper").eq(0);
                wrapper.addClass("current");
                $(".help-dropdown", wrapper).slideToggle();
                var otherHelpWrapper = $(".help-wrapper:not(.current)");
                $(".help-dropdown", otherHelpWrapper).slideUp();
                wrapper.removeClass("current");
            });
            //hide help dropdown if click ouside 
            $("body").click(function (e) {
                if ($(".help-dropdown:visible").length > 0) {
                    if ($(e.target).closest('.help-wrapper').length === 0) {
                        $(".help-dropdown:visible").slideUp();
                    }
                }
            });
            //file uploader
            $(".file-uploader-wrapper a").click(function () {
                var wrapper = $(this).parents(".file-uploader-wrapper").eq(0);
                $('input.text-file', wrapper).eq(0).trigger('click');
            });

            //clear selected file
            $(".file-uploader-wrapper input.text").keypress(function (e) {
                var charCode = (e.which) ? e.which : e.keyCode;
                if (charCode === 8) {
                    var wrapper = $(this).parents(".file-uploader-wrapper").eq(0);
                    $('input.text-file', wrapper).eq(0).val("").trigger('change');
                    e.preventDefault();
                }
            });

            //radio input click
            $(".radio-wrapper label").click(function () {
                var wrapper = $(this).parents(".radio-wrapper").eq(0);
                $(".input-dot", wrapper).removeClass("input-dot-checked");
                var forId = $(this).attr("for");
                var input = $("input#" + forId);

                if (input.hasClass("radio")) {
                    $(".input-dot[for='" + forId + "']", wrapper).addClass("input-dot-checked").find("a").focus();
                } else if (input.hasClass("checkbox")) {
                    if (input.prop("checked")) {
                        $(".input-dot[for='" + forId + "']", wrapper).addClass("input-dot-checked").find("a").focus();
                    } else {
                        $(".input-dot[for='" + forId + "']", wrapper).find("a").focus();
                    }
                }

                if (input.hasClass("js-disable-fieldset")) {
                    var tId = input.attr('data-ctrl-target');
                    var wrapper = $("." + tId);
                    $("input.text", wrapper).val("").prop("disabled", true).removeClass("input-error"); ;
                    $("fieldset", wrapper).addClass("fieldset-disabled");

                } else if (input.hasClass("js-enable-fieldset")) {
                    var tId = input.attr('data-ctrl-target');
                    var wrapper = $("." + tId);
                    $("input.text", wrapper).val("").prop("disabled", false).removeClass("input-error"); ;
                    $("fieldset", wrapper).removeClass("fieldset-disabled");
                }
                else if (input.hasClass("js-toggle-fieldset")) {
                    if (input.prop("checked")) {
                        var tId = input.attr('data-ctrl-target');
                        var wrapper = $("." + tId);
                        $("input.text", wrapper).each(function () {
                            $(this).prop("disabled", false).removeClass("input-error");
                            if ($(this)[0].hasAttribute("data-default-value")) {
                                $(this).val($(this).attr("data-default-value"));
                            } else {
                                $(this).val("");
                            }
                        });
                        wrapper.removeClass("fieldset-disabled");
                    } else {
                        var tId = input.attr('data-ctrl-target');
                        var wrapper = $("." + tId);
                        $("input.text", wrapper).val("").prop("disabled", true).removeClass("input-error"); ;
                        wrapper.addClass("fieldset-disabled");

                    }
                }

                window.setTimeout(function () {
                    $(".input-dot[for='" + forId + "'] a", wrapper).focus();
                }, 200);
            });

            //trigger radio/checkbox label click
            $(".radio-wrapper .input-dot a").click(function () {
                $(this).parent().trigger("click");
            });

            //keypress a radio dot, trigger click
            $(".radio-wrapper .iput-dot a").keypress(function (e) {
                if (e.keyCode === 13) { // If the the enter key was pressed.
                    $(this).trigger("click");
                }
            });

            //click unit to focus input box
            $(".input-wrapper .unit").click(function () {
                $(this).prev("input.text").focus();
            });

            //focus a text input
            $(".input-wrapper input.text").keyup(function () {
                if ($(this).val() === "") {
                    $(this).next(".unit-left").show();
                } else {
                    $(this).next(".unit-left").hide();
                }
            });

            //allow only valid char for filename
            $(".text-filename, .text-filename-without-extension").keypress(function (e) {
                return openWARP.app.isValidCharForFileName(e);
            });

            //validate filename if paster in invalid char!!
            $(".text-filename, .text-filename-without-extension").keyup(function (e) {
                var hasExtension = !$(this).hasClass("text-filename-without-extension");
                if ($(this).val().length === 0 || openWARP.app.isFileName($(this).val(), hasExtension)) {
                    $(this).removeClass("input-error");
                } else {
                    $(this).addClass("input-error");
                }
            });
            $(".text-filename, .text-filename-without-extension").blur(function (e) {
                if ($(this).hasClass("text-filename-without-extension")) {
                    $(this).trigger("keyup");
                } else {
                    if ($(this).val().length === 0 || openWARP.app.isFileNameWithExtension($(this).val())) {
                        $(this).removeClass("input-error");
                    } else {
                        $(this).addClass("input-error");
                    }
                }
            });
            $(".text-filename, .text-filename-without-extension").focus(function (e) {
                $(this).trigger("keyup");
            });

            //allow only numic been input
            $("input.text-num").keypress(function (e) {
                return openWARP.app.isValidCharForNum(e, $(this).val());
            });

            //validate numric format correct
            $("input.text-num").keyup(function (e) {
                var isAllowE = $(this).hasClass("text-e-num");
                if ($(this).val().length === 0 || openWARP.app.isNumic($(this).val(), isAllowE) || $(this).hasClass("text-list-num")) {
                    $(this).removeClass("input-error");
                } else {
                    $(this).addClass("input-error");
                }
                if ($(this).val().length > 0) {
                    var max = $(this).attr("data-max-value");
                    var min = $(this).attr("data-min-value");
                    var cVal = parseFloat($(this).val());
                    if ((typeof max !== typeof undefined) && (max !== false)) {
                        max = parseFloat(max);
                        if (cVal > max) {
                            $(this).addClass("input-error");
                        }
                    }
                    if ((typeof min !== typeof undefined) && (min !== false)) {
                        min = parseFloat(min);
                        if (cVal < min) {
                            $(this).addClass("input-error");
                        }
                    }
                }
                if ($(this).val() === "-") {
                    $(this).removeClass("input-error");
                }
            });
            $("input.text-num").blur(function (e) {
                $(this).trigger("keyup");
            });
            $("input.text-num").focus(function (e) {
                $(this).trigger("keyup");
            });

            //validate input value after data loaded from local storage
            window.setTimeout(function () {
                $(" input.text-filename-without-extension:visible, input.text-num:visible").trigger("keyup");
                $("input.text-filename:visible").trigger("blur");
                $(".input-wrapper .unit-left:visible").each(function () {
                    if ($(this).prev("input.text").val() === "") {
                        $(this).show();
                    } else {
                        $(this).hide();
                    }
                });
            }, 100);
        },

        // Used to construct 3-dimension vector from HTML elements.
        construct3dVectorParameter: function (isTranslate, keyName, axis) {
            var line = "";
            line += isTranslate ? "1 " : "2 ";
            if (axis == 0) {
                line += "1. 0. 0. ";
            } else if (axis == 1) {
                line += "0. 1. 0. ";
            } else {
                line += "0. 0. 1. ";
            }
            var invalidate = false;
            $.each([1, 2, 3], function (i, item) {
                var item_val = $(keyName + "_" + item).val();
                if (!item_val || item_val == "") {
                    invalidate = true;
                    return false;
                }
                line += (item_val + " ");
            });
            if (invalidate) {
                return "";
            }
            return line;
        },

        // Used to pop up a dialog to show specified message.
        message_box: function (label, content) {
            // Remove the loading features
            openWARP.app.hide_loading();
            $("#dialog_label").html(label);
            $("#dialog_content").val(content);
            $("#single_dialog").dialog("open");
        },

        // Used to pop up a loading dialog to tell user that the application is processing the request
        show_loading: function () {
            $("#dialog_label").html("Loading ...");
            $("#dialog_content").hide();
            $(".ui-dialog-titlebar-close", $("#single_dialog").parent()).hide();
            $("#single_dialog").dialog("open");
        },

        // Used to hide the loading dialog
        hide_loading: function () {
            $("#dialog_content").show();
            $(".ui-dialog-titlebar-close", $("#single_dialog").parent()).show();
            $("#single_dialog").dialog("close");
        },

        // Used to quit the application
        quit_app: function () {
            $.ajax({
                url: "/quit",
                type: "POST",
                success: function (data) {
                    // Instead of closing browser, we redirect the page to a goodbye page.
                    window.location = "/bye.html";
                },
                error: function (request, status, error) {
                    // Notify user
                    openWARP.app.message_box("Error in server side:", error + "\n\nDetails:\n" + request.responseJSON.error);
                }
            });
        },

        isValidCharForFileName: function (evt) {
            var charCode = (evt.which) ? evt.which : evt.keyCode
            var str = '\\/:*?"<>|'; // forbidden characters \ / : * ? " < > | 
            if ($(evt.target).hasClass("text-filename-without-extension")) {
                str += "."; // forbidden characters . if no extension
            }
            for (var i = 0; i < str.length; i++) {
                var n = str.charCodeAt(i);
                if (charCode === n) {
                    return false;
                }
            }
            return true;
        },
        isFileName: function (str, hasExtension) {
            if (hasExtension) {
                var regex = /^[^\\\/:\*\?"<>\|]+$/; // forbidden characters \ / : * ? " < > | .
            } else {
                var regex = /^[^\\\/:\*\?"<>\|.]+$/; // forbidden characters \ / : * ? " < > | .
            }
            return regex.test(str);
        },
        isFileNameWithExtension: function (str) {
            var regex = /^([^\\\/:\*\?"<>\|])+(\.)+([a-zA-Z0-9]{2,4})$/;
            return regex.test(str);

            //var regex = /^([a-zA-Z0-9_.+-])+\@(([a-zA-Z0-9-])+\.)+([a-zA-Z0-9]{2,4})+$/;

        },
        isValidCharForNum: function (evt, cVal) {
            var ctrlDown = false;
            evt = evt || window.event // IE support
            var ctrlDown = evt.ctrlKey || evt.metaKey // Mac support
            // Check for Alt+Ctrl
            if (ctrlDown && evt.altKey) return false;
            var charCode = (evt.which) ? evt.which : evt.keyCode
            //allows CTRL + c/v/x
            //c			99
            //C			67
            //v			118
            //V			86
            //x			120
            //X			88
            //left arrow		37 
            //right arrow		39 
            //allows char . -
            var str = ".-eE";
            var n1 = str.charCodeAt(0);
            var n2 = str.charCodeAt(1);
            var n3 = str.charCodeAt(2);
            var n4 = str.charCodeAt(3);

            var isAllowE = false;
            if ($(evt.target).hasClass("text-e-num")) {
                var isAllowE = (charCode === n3) || (charCode === n4);
            }

            if ($(evt.target).hasClass("text-list-num")) {
                var isAllowE = (charCode === ' '.charCodeAt(0));
            }

            if (charCode > 31 && (charCode < 48 || charCode > 57) && (charCode !== 37) && (charCode !== 39) && (charCode !== n1) && (charCode !== n2) && (!isAllowE)) {
                if (!(ctrlDown && charCode === 99) && !(ctrlDown && charCode === 67) && !(ctrlDown && charCode === 118) && !(ctrlDown && charCode === 86) && !(ctrlDown && charCode === 120) && !(ctrlDown && charCode === 88))
                    return false;
            }

            //check for two dots
            cVal = String(cVal);
            if (charCode === n1 && (cVal.indexOf(".") >= 0)) {
                return false;
            }
            //check for two minus
            if (charCode === n2 && (cVal.indexOf("-") >= 0)) {
                return false;
            }
            //add zero if directly input .
            if (charCode === n1 && (cVal.length === 0)) {
                $(evt.target).val("0")
            }
            //add zero if directly input -.
            if (charCode === n1 && (cVal === "-")) {
                $(evt.target).val("-0")
            }
            return true;
        },
        isNumic: function (str, isAllowE) {
            var isNum = false;
            if ($.isNumeric(str) && (str.indexOf("x") < 0) && (str.indexOf("X") < 0)) {
                isNum = true;
            }
            if (!isAllowE && ((str.indexOf("e") >= 0) || (str.indexOf("E") >= 0))) {
                isNum = false;
            }
            if ((str.length > 1) && (str.indexOf("0") === 0) && (str.indexOf(".") !== 1)) {
                isNum = false;
            }
            return isNum;
        },
        // Check if test_str is a non-negative integer.
        isNonNegativeInteger: function (test_str) {
            var num = parseInt(test_str);
            if ("" + num != test_str || num < 0) {
                return false;
            }
            return true;
        }
    }
}

//Display filename after select file.
function uploadFileTrigger(obj) {
    // Calculate the file name
	var fileUploaderWrapper = $(obj).parent();
	var url = $(obj).val();
	var lastIndex = url.lastIndexOf('\\');
	var fileName = url.substring(lastIndex + 1);
	var lastIndex = url.lastIndexOf('/');
	var name2 = url.substring(lastIndex + 1);
	if (fileName.length > name2.length) {
		fileName = name2;
	}
    // Set file name to the text element
	$(obj).parent().find("input.text").val(fileName).trigger("change");
	var obj_id = $(obj).attr("id");
	// Upload the file
	if (url.length > 0) {
	    $.ajaxFileUpload({
	        url: "upload_file",
	        secureuri: false,
	        fileElementId: $(obj).attr("id"),
	        dataType: "text/json",
	        success: function (data, status) {
	            // Update the hidden filepath field
                var json_data = jQuery.parseJSON(jQuery(data).text());
                var filepath = json_data.filepath;
	            var wrapper = $("#" + obj_id).parents(".file-uploader-wrapper");
	            $(":hidden", wrapper).val(filepath);

	            // Update some data of the page if necessary
	            if ($("#" + obj_id).hasClass("js-ctrl-fieldset")) {
	                var tId = $("#" + obj_id).attr('data-ctrl-target');
	                var wrapper_points_panels = $("." + tId);
	                if (filepath.length < 1 || !json_data.points || !json_data.panels) {
	                    $("input.text", wrapper_points_panels).val("").prop("disabled", true).removeClass("input-error");
	                    $(".unit", wrapper_points_panels).show();
	                    wrapper_points_panels.addClass("fieldset-disabled");
	                } else {
	                    // Points&panels value come from uploaded mesh file analysis
	                    $("input.text", wrapper_points_panels).eq(0).val(json_data.points).prop("disabled", false).removeClass("input-error");
	                    $("input.text", wrapper_points_panels).eq(1).val(json_data.panels).prop("disabled", false).removeClass("input-error");
	                    $(".unit", wrapper_points_panels).hide();
	                    wrapper_points_panels.removeClass("fieldset-disabled");
	                }
	            }
	        },
	        error: function (data, status, error) {
	            openWARP.app.message_box("File upload error:", ""+error);
	        }
	    });
    }
}

/*
*  document ready
*  This function will render when the document is ready (page has loaded)
*/
 
$(document).ready(function(){
	//call functions for all pages 
     openWARP.common.init();

	 //call functions for landing page
	if ($(".landing-wrapper").length > 0){
		openWARP.landing.init();
	}

	//call functions for app page
	if ($(".openwarp-app").length > 0){
		openWARP.app.init();
	}
});