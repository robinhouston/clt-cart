#!/usr/bin/env node

var d3 = require("d3"),
    fs = require("fs"),
    topojson = require("topojson");

function write_csv_row(row) {
	console.log(row.map(function(el) {
		if (("" + el).match(/[",\n]/)) {
			return '"' + el.replace(/"/g, '""') + '"';
		} else {
			return el;
		}
	}).join(","));
}

var path = d3.geo.path().projection(null);

if (process.argv.length == 2) {
	console.error("No topojson file supplied");
	process.exit(1);
}
else if (process.argv.length > 3) {
	console.error("Unexpected argument: " + process.argv[3]);
	process.exit(1);
}
var topojson_filename = process.argv[2];
var t = JSON.parse(fs.readFileSync(topojson_filename))

for (var object_name in t.objects) {
	var feature_collection = topojson.feature(t, t.objects[object_name]),
	    features = feature_collection.features;
	write_csv_row(["Object", "Feature ID", "Area"]);
	features.forEach(function(feature) {
		write_csv_row([object_name, feature.id, path.area(feature)]);
	});
}
