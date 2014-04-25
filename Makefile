out/cart/areas.txt: out/cart/map.topo.json
	bin/topojson-area out/cart/map.topo.json > $@

out/cart/map.topo.json: data/map.topo.json out/cart/points.txt
	bin/tj-update-points.py data/map.topo.json out/cart/points.txt > $@

out/cart/points.txt: data/density.grid data/points.txt
	code/cart 500 250 data/density.grid data/points.txt > $@

data/density.grid: data/population.csv data/area.csv
	bin/density-grid.py data/map-grid.csv \
	    --numerator-key=Country --numerator-value=Population data/population.csv \
	    --denominator-key="Feature ID" --denominator-value=Area data/area.csv \
	> $@

data/area.csv: data/map.topo.json
	bin/topojson-area data/map.topo.json > $@

data/points.txt: data/map.topo.json
	bin/tj-extract-points.py data/map.topo.json > $@

data/map.topo.json: data/map.geo.json
	topojson -o $@ data/map.geo.json

data/map.geo.json:
	simplification=5000
	/usr/local/cartograms/bin/as-js.py --format=geojson --map=world-robinson \
		--simplification="$simplification" > $@

data/population.csv:
	(
	    echo "Country,Population"
	    /usr/local/cartograms/bin/csv-snip --cols=B,E --expect-header="Alpha-2,Value" \
	    --munge=E=s/\.0$// \
	    ~/Kiln/world-gold-council/GoldMap/data/maps/live/Population.csv  | sort
	) > $@
