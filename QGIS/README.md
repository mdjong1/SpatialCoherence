# How to make spatial coherence map in QGIS

- Load TIF from SpatialCoherence code
- Raster -> Conversion -> Polygonize (to create polygons from each raster cell, band irrelevant)
- Vector -> Geometry -> Centroids (on created vectorized layer)
- Point sample these centroids on original TIF using Point Sampling plugin
- Sample for each band (1 = entry, 2 = exit, 3 = active) and save
- Copy and paste entry and exit style for respective layers
- Adjust color range depending on min/max values of these results in layer Symbology -> Simple Marker -> Fill color EQUATION -> Change large number on scale_linear
- Order layers correctly to have smaller squares on top (Exits above Entries) or switch styles
