---
layout: page
title: About
---

<p class="message">
  Hey there! This page is included as an example. Feel free to customize it for your own use upon downloading. Carry on!
</p>

In the novel, *The Strange Case of Dr. Jeykll and Mr. Hyde*, Mr. Poole is Dr. Jekyll's virtuous and loyal butler. Similarly, Poole is an upstanding and effective butler that helps you build Jekyll themes. It's made by [@mdo](https://twitter.com/mdo).

There are currently two themes built on Poole:

* [Hyde](http://hyde.getpoole.com)
* [Lanyon](http://lanyon.getpoole.com)

Learn more and contribute on [GitHub](https://github.com/poole).

## Setup

Some fun facts about the setup of this project include:

* Built for [Jekyll](http://jekyllrb.com)
* Developed on GitHub and hosted for free on [GitHub Pages](https://pages.github.com)
* Coded with [Sublime Text 2](http://sublimetext.com), an amazing code editor
* Designed and developed while listening to music like [Blood Bros Trilogy](https://soundcloud.com/maddecent/sets/blood-bros-series)

Have questions or suggestions? Feel free to [open an issue on GitHub](https://github.com/poole/issues/new) or [ask me on Twitter](https://twitter.com/mdo).

Thanks for reading!

<head>
    <title>Simple Leaflet Map</title>
    <meta charset="utf-8" />
    <link 
        rel="stylesheet" 
        href="http://cdn.leafletjs.com/leaflet-0.7/leaflet.css"
    />
</head>
<body>
    <div id="map" style="width: 600px; height: 400px"></div>

    <script
        src="http://cdn.leafletjs.com/leaflet-0.7/leaflet.js">
    </script>

    <script>
        var map = L.map('map').setView([-41.2858, 174.78682], 14);
        mapLink = 
            '<a href="http://openstreetmap.org">OpenStreetMap</a>';
        L.tileLayer(
            'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; ' + mapLink + ' Contributors',
            maxZoom: 18,
            }).addTo(map);
    </script>
</body>
