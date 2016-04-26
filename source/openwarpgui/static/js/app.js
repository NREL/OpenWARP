'use strict';
/*jshint -W097*/
/*jshint strict:false*/
/*global angular, console, document*/

// main angular configurations
var openwarp = angular.module('openwarp', ['openwarp.controllers', 'ngRoute', 'LocalStorageModule']);

// configure routes
openwarp.config(['$routeProvider', function ($routeProvider) {
	$routeProvider
        .when('/meshing', {
            templateUrl: 'partials/meshing.html',
            controller: 'meshing'
        })
        .when('/simulation', {
            templateUrl: 'partials/simulation.html',
            controller: 'simulation'
        })
        .when('/postprocessing', {
            templateUrl: 'partials/postprocessing.html',
            controller: 'postprocessing'
        })
        .when('/configuration', {
            templateUrl: 'partials/logging.html',
            controller: 'configuration'
        })
        .otherwise({
            redirectTo: '/meshing'
        });
}]).config(['localStorageServiceProvider', function(localStorageServiceProvider){
  localStorageServiceProvider.setPrefix('openWARP');
}]);