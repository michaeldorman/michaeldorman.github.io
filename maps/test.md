
---
layout: default
title: About
---

<div id="map">

</div>

<script>
        var map = L.map('map').setView([34.00000, -118.260126], 14); 

        mapLink =
'<a href="http://openstreetmap.org">OpenStreetMap</a>'; L.tileLayer(
'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', { attribution: 'Map data &copy; ' + mapLink, maxZoom: 18, }).addTo(map);


        var marker = L.marker([34.063298, -118.260126]) .addTo(map).bindPopup("<b>Blah Blah Blah</b><br /><a href='http://www.cnn.com'>Additional Information</a><br />").openPopup();


</script>
