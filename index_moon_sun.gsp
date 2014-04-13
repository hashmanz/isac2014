

<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>Gravity Map</title>
    <style>
      html, body, #map-canvas {
        height: 100%;
        margin: 0px;
        padding: 0px
      }
      #panel {
        position: absolute;
        top: 5px;
        left: 50%;
        margin-left: -180px;
        z-index: 5;
        background-color: #fff;
        padding: 5px;
        border: 1px solid #999;
      }
      section {
	    width: 100%;
	    height: 100%;
	    margin: auto;
	    padding: 0px;
	  }
		div#one {
		    width: 15%;
		    height: 200px;
		    float: left;
		}
		div#map-canvas {
		    margin-left: 15%;
		    height: 100%;
		    background: #2E64FE;
		}
    </style>
    <style type="text/css">
	   .labels {
	     color: black;
	     font-family: "Lucida Grande", "Arial", sans-serif;
	     font-size: 10px;
	     text-align: center;
	     width: 10px;     
	     white-space: nowrap;
	   }
	 </style>
	 <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js" type="text/javascript"></script>
     <script src="https://maps.googleapis.com/maps/api/js?v=3.exp&sensor=false&libraries=visualization"></script>
     <script type="text/javascript">
		var map, pointarray, heatmap;
		
		
		var taxiData = [
		  {location: new google.maps.LatLng(-34.782551, 151.445368), weight: 10},
		  {location: new google.maps.LatLng(-34.782551, 155), weight: 10}
		  //new google.maps.LatLng(-34.782551, 151.445368),
		  //new google.maps.LatLng(-34.782551, 151.445368),
		  //new google.maps.LatLng(-34.782551, 151.445368)
		
		];

		var gravityData=[];
    	var gravities_around = "${gravities_around}".replace("[","").replace("]","").split(",");
    	var lats_around = "${lats_around}".replace("[","").replace("]","").split(",");
    	var lons_around = "${lons_around}".replace("[","").replace("]","").split(",");
    
    	for (var i=0; i<gravities_around.length; i++) {
        	lat = parseFloat(lats_around[i]);
        	lon = parseFloat(lons_around[i]);
        	w = parseFloat(gravities_around[i]);
    		gravityData[i] = {location: new google.maps.LatLng(lat, lon), weight: w};
    		//gravityData[i] = new google.maps.LatLng(lats_around[i], lons_around[i]);
    	}

		
    	
		$( document ).ready(function() {
		//function initialize(){
			var lat = ${lat};
			var lon = ${lon};
			
			var userLatLng = new google.maps.LatLng(lat,lon );

			// Create ElevationService
			elevator = new google.maps.ElevationService();
			var obj=new Object();
			obj.latLng=userLatLng;
			getElevation(obj);
			
//			if(!navigator.geolocation) {
//				navigator.geolocation.getCurrentPosition(function(position) {
//					//position = navigator.geolocation.getCurrentPosition();
//					userLatLng = new google.maps.LatLng(position.coords.latitude, position.coords.longitude);
//					lat = position.coords.latitude;
//					lon = position.coords.longitude
//				});
//			}

			if (lat>=0){
				$("#lat_direction").val("North");
				$("#lat_value").val(lat);
			}
			if (lat<0){
				$("#lat_direction").val("South");
				$("#lat_value").val(-lat);
			}
			if (lon>=0){
				$("#lon_direction").val("East");
				$("#lon_value").val(lon);
			}
			if (lon<0){
				$("#lon_direction").val("West");
				$("#lon_value").val(-lon);
			}
					
			var mapOptions = {
			  zoom: 12,
			  center: userLatLng,
			  mapTypeId: google.maps.MapTypeId.SATELLITE
			  //mapTypeId: google.maps.MapTypeId.ROADMAP
			};
		
			map = new google.maps.Map(document.getElementById('map-canvas'), mapOptions);
		      
		 
			var marker = new google.maps.Marker({
			    position: userLatLng,
			    map: map,
			    draggable: true,
			       raiseOnDrag: true,
			       labelClass: "labels", // the CSS class for the label
			       labelInBackground: false,
			       labelStyle: {opacity: 0.75}
			});
				
			var gravity = "${gravity}";
			var iw = new google.maps.InfoWindow({
			    content: "Gravity=${gravity}"
			});
			
			google.maps.event.addListener(marker, "click", function (e) { iw.open(map, this); });
			iw.open(map,marker);


			google.maps.event.addListenerOnce(map, 'idle', function(){
				var bounds = this.getBounds();
		        var ne = bounds.getNorthEast(); // LatLng of the north-east corner
				var sw = bounds.getSouthWest(); // LatLng of the south-west corder
				var nw = new google.maps.LatLng(ne.lat(), sw.lng());
				var se = new google.maps.LatLng(sw.lat(), ne.lng());
				$('input[name="ne_lat"]').val(ne.lat());
				$('input[name="ne_lon"]').val(ne.lng());
				$('input[name="sw_lat"]').val(sw.lat());
				$('input[name="sw_lon"]').val(sw.lng());
				$('input[name="nw_lat"]').val(nw.lat());
				$('input[name="nw_lon"]').val(nw.lng());
				$('input[name="se_lat"]').val(se.lat());
				$('input[name="se_lon"]').val(se.lng());
		    });

			google.maps.event.addListener(map, 'bounds_changed', updateBoundarys);

			
			
			var pointArray = new google.maps.MVCArray(gravityData);
			
			heatmap = new google.maps.visualization.HeatmapLayer({
			  data: pointArray
			});
			
			heatmap.setMap(map);


			google.maps.event.addListener(map, 'click', function(event) {
				//if (drag){return;}
				//posset = 1;
				//fc(event.latLng) ;
				//if (map.getZoom() < 10){map.setZoom(10);}
				map.panTo(event.latLng);
				computepos(event.latLng);
			});
		});
		//}

		function getElevation(event) 
		{
			var locations = [];
	
			// Retrieve the clicked location and push it on the array
			var clickedLocation = event.latLng;
			locations.push(clickedLocation);
	
			// Create a LocationElevationRequest object using the array's one value
			var positionalRequest = {'locations': locations};
	
			// Initiate the location request
			elevator.getElevationForLocations(positionalRequest, function(results, status) 
			{
				if (status == google.maps.ElevationStatus.OK) 
				{
					// Retrieve the first result
					if (results[0]) 
					{
						alt = Math.round(results[0].elevation)
						$('input[name="alt_value"]').val(alt);
					}
				}
			});
		}

		function computepos (point)
		{
			var latA = Math.abs(Math.round(point.lat()));
			
			var lonA = Math.abs(Math.round(point.lng()));
			

			if(point.lat() < 0)
			{
				$("#lat_direction").val("South");
				$("#lat_value").val(latA);
				$('input[name="lat_value"]').val(-latA);
			}
			else
			{
				$("#lat_value").val(latA);
				$('input[name="lat_value"]').val(latA);
			}
			if(point.lng() < 0)
			{
				$("#lon_direction").val("West");
				$("#lon_value").val(lonA);
				$('input[name="lon_value"]').val(-lonA);
			}
			else
			{
				$("#lon_value").val(lonA);
				$('input[name="lon_value"]').val(lonA);
			}
		}

		function updateBoundarys(event){
			var bounds = this.getBounds();
	        var ne = bounds.getNorthEast(); // LatLng of the north-east corner
			var sw = bounds.getSouthWest(); // LatLng of the south-west corder
			var nw = new google.maps.LatLng(ne.lat(), sw.lng());
			var se = new google.maps.LatLng(sw.lat(), ne.lng());
			$('input[name="ne_lat"]').val(ne.lat());
			$('input[name="ne_lon"]').val(ne.lng());
			$('input[name="sw_lat"]').val(sw.lat());
			$('input[name="sw_lon"]').val(sw.lng());
			$('input[name="nw_lat"]').val(nw.lat());
			$('input[name="nw_lon"]').val(nw.lng());
			$('input[name="se_lat"]').val(se.lat());
			$('input[name="se_lon"]').val(se.lng());
	    }
		
		function toggleHeatmap() {
			heatmap.setMap(heatmap.getMap() ? null : map);
		}
		
		function changeGradient() {
			var gradient = [
			  'rgba(0, 255, 255, 0)',
			  'rgba(0, 255, 255, 1)',
			  'rgba(0, 191, 255, 1)',
			  'rgba(0, 127, 255, 1)',
			  'rgba(0, 63, 255, 1)',
			  'rgba(0, 0, 255, 1)',
			  'rgba(0, 0, 223, 1)',
			  'rgba(0, 0, 191, 1)',
			  'rgba(0, 0, 159, 1)',
			  'rgba(0, 0, 127, 1)',
			  'rgba(63, 0, 91, 1)',
			  'rgba(127, 0, 63, 1)',
			  'rgba(191, 0, 31, 1)',
			  'rgba(255, 0, 0, 1)'
			]
			heatmap.set('gradient', heatmap.get('gradient') ? null : gradient);
		}
		
		function changeRadius() {
		  	heatmap.set('radius', heatmap.get('radius') ? null : 20);
		}
			
		function changeOpacity() {
		  	heatmap.set('opacity', heatmap.get('opacity') ? null : 0.2);
		}
			
		google.maps.event.addDomListener(window, 'load', initialize);


	    $( document ).ready(function() {

	    	$("#lat_direction").change(function () {
	            var lat_direction = this.value;
	            var lat_value = $("#lat_value").val();
	            if (lat_direction=="South"){
	            	lat_value = -lat_value;
		        }
	            $('input[name="lat_value"]').val(lat_value);
		    });
	    	$("#lat_value").change(function () {
	            var lat_value = this.value;
	            var lat_direction = $("#lat_direction").val();
	            if (lat_direction=="South"){
	            	lat_value = -lat_value;
		        }
	            $('input[name="lat_value"]').val(lat_value);
		    });
	    	$("#lon_direction").change(function () {
	            var lon_direction = this.value;
	            var lon_value = $("#lon_value").val();
	            if (lon_direction=="West"){
	            	lon_value = -lon_value;
		        }
	            $('input[name="lon_value"]').val(lon_value);
		    });
	    	$("#lon_value").change(function () {
	            var lon_value = this.value;
	            var lon_direction = $("#lon_direction").val();
	            if (lon_direction=="South"){
	            	lon_value = -lon_value;
		        }
	            $('input[name="lon_value"]').val(lon_value);
		    });

	    	
	    });
		
// Leon's Moon/Sun Code
// This code was developed to calculate the Moon and Sun Gravity at any Point
// This involved modifying and developing existing code and formulas from the following links .
// Reference: http://www.dept.aoe.vt.edu/~lutze/AOE4134/13LocalSiderealTime.pdf
// http://www.bogan.ca/astro/telescopes/coodcvtn.html
// http://www.lunar-occultations.com/rlo/ephemeris.htm - developed by Keith Burnett

   	function doCalcs(form) {
   	var g, days,t ,L1, M1, C1, V1, Ec1, R1, Th1, Om1, Lam1, Obl, Ra1, Dec1;
	var F, L2, Om2, M2, D, R2, R3, Bm, Lm, HLm, HBm, Ra2, Dec2, EL, EB, W, X, Y, A; 
	var Co, SLt, Psi, Il, K, P1, P2, y, m, d, bit, h, min, bk;
	var Latitude, Longitude, LatRaw, LongRaw, MajorAxis, MinorAxis, WGSHeight, XEllipse, YEllipse, RadiusEarth, rightNow
	var JD, Tut1, UT, LocSidTime, UTSidModTime, UTSidTime, LocHourAngleMoon, LocHourAngleSun, SineAltitudeMoon, SineAltitudeSun, JuDay
	var EARTHTOMOON, MASSMOON, GCONST, aQuad, bQuad, cQuad, DistMoon, EARTHTOSUN, MASSSUN
	var GrossForceMoon, NormalForceMoon, GrossForceSun, NormalForceSun
//
//	Get date and time code from user, isolate the year, month, day and hours
//	and minutes, and do some basic error checking! This only works for AD years
//  
//	g = form.num.value;
//	y = Math.floor(g / 10000);
//	m = Math.floor( (g - y*10000) / 100);
//	d = Math.floor(g - y*10000 - m*100 );
//	bit = (g - Math.floor(g))*100;
//	h = Math.floor(bit);
//	min = Math.floor(bit*100 - h * 100 + 0.5);
	

	rightNow = new Date();
	y = rightNow.getUTCFullYear();
	m = rightNow.getUTCMonth() + 1;
	d = rightNow.getUTCDate();
	h = rightNow.getUTCHours();
	min = rightNow.getUTCMinutes();
	
//
//	primative error checking - accounting for right number of
//	days per month including leap years. Using bk variable to 
//	prevent multiple alerts. See functions isleap(y) 
//	and goodmonthday(y, m, d).
//
	bk = 0;
	if(g < 16000000) {	
		bk = 1;
		alert("Routines are not accurate enough to work back that" + 
               " far - answers are meaningless!");
		}
	if(g > 23000000) {
		bk = 1;	
		alert("Routines are not accurate enough to work far into the future" + 
               " - answers are meaningless!");
		}	
	if( ((m<1) || (m>12)) && (bk !=1)) {
			bk = 1;
			alert("Months are not right - type date again");
			}
	if((goodmonthday(y, m, d) ==0) && (bk !=1)) {
			bk = 1;
			alert("Wrong number of days for the month or not a leap year - type date again");
			}
	if ((h>23) && (bk != 1)) {
			bk = 1;
			alert("Hours are not right - type date again");
			}
	if((min>59) && (bk !=1)) {
			alert("Minutes are not right - type date again");
			}
//
//	Get the number of days since J2000.0 using day2000() function
//
    days = day2000(y, m, d, h + min/60);
	t = days / 36525;

//
//	Sun formulas
//
//	L1	- Mean longitude
//	M1	- Mean anomaly
//	C1	- Equation of centre
//	V1	- True anomaly
//	Ec1	- Eccentricity 
//	R1	- Sun distance
//	Th1	- Theta (true longitude)
//	Om1	- Long Asc Node (Omega)
//	Lam1- Lambda (apparent longitude)
//	Obl	- Obliquity of ecliptic
//	Ra1	- Right Ascension
//	Dec1- Declination
//

	L1 = range(280.466 + 36000.8 * t);
	M1 = range(357.529+35999*t - 0.0001536* t*t + t*t*t/24490000);
	C1 = (1.915 - 0.004817* t - 0.000014* t * t)* dsin(M1);	 
	C1 = C1 + (0.01999 - 0.000101 * t)* dsin(2*M1);
	C1 = C1 + 0.00029 * dsin(3*M1);
	V1 = M1 + C1;
	Ec1 = 0.01671 - 0.00004204 * t - 0.0000001236 * t*t;
	R1 = 0.99972 / (1 + Ec1 * dcos(V1));
	Th1 = L1 + C1;
	Om1 = range(125.04 - 1934.1 * t);
	Lam1 = Th1 - 0.00569 - 0.00478 * dsin(Om1);
	Obl = (84381.448 - 46.815 * t)/3600;
	Ra1 = datan2(dsin(Th1) * dcos(Obl) - dtan(0)* dsin(Obl), dcos(Th1));
	Dec1 = dasin(dsin(0)* dcos(Obl) + dcos(0)*dsin(Obl)*dsin(Th1));

//
//	Moon formulas
//
//	F 	- Argument of latitude (F)
//	L2 	- Mean longitude (L')
//	Om2 - Long. Asc. Node (Om')
//	M2	- Mean anomaly (M')
//	D	- Mean elongation (D)
//	D2	- 2 * D
//	R2	- Lunar distance (Earth - Moon distance)
//	R3	- Distance ratio (Sun / Moon)
//	Bm	- Geocentric Latitude of Moon
//	Lm	- Geocentric Longitude of Moon
//	HLm	- Heliocentric longitude
//	HBm	- Heliocentric latitude
//	Ra2	- Lunar Right Ascension
//	Dec2- Declination
//

	F = range(93.2721 + 483202 * t - 0.003403 * t* t - t * t * t/3526000);
	L2 = range(218.316 + 481268 * t);
	Om2 = range(125.045 - 1934.14 * t + 0.002071 * t * t + t * t * t/450000);
	M2 = range(134.963 + 477199 * t + 0.008997 * t * t + t * t * t/69700);
	D = range(297.85 + 445267 * t - 0.00163 * t * t + t * t * t/545900);
	D2 = 2*D;
	R2 = 1 + (-20954 * dcos(M2) - 3699 * dcos(D2 - M2) - 2956 * dcos(D2)) / 385000;
	R3 = (R2 / R1) / 379.168831168831;
	Bm = 5.128 * dsin(F) + 0.2806 * dsin(M2 + F);
	Bm = Bm + 0.2777 * dsin(M2 - F) + 0.1732 * dsin(D2 - F);
	Lm = 6.289 * dsin(M2) + 1.274 * dsin(D2 -M2) + 0.6583 * dsin(D2); 
	Lm = Lm + 0.2136 * dsin(2*M2) - 0.1851 * dsin(M1) - 0.1143 * dsin(2 * F); 
	Lm = Lm +0.0588 * dsin(D2 - 2*M2) 
	Lm = Lm + 0.0572* dsin(D2 - M1 - M2) + 0.0533* dsin(D2 + M2);
	Lm = Lm + L2;
	Ra2 = datan2(dsin(Lm) * dcos(Obl) - dtan(Bm)* dsin(Obl), dcos(Lm));
	Dec2 = dasin(dsin(Bm)* dcos(Obl) + dcos(Bm)*dsin(Obl)*dsin(Lm));
	HLm = range(Lam1 + 180 + (180/Math.PI) * R3 * dcos(Bm) * dsin(Lam1 - Lm));
	HBm = R3 * Bm;

// INPUTS: LAT, LONG, ALT
// WGS84 major and minor axes
// Assume height is WGS84 height
    LatRaw = $('input[name="lat_value"]').val();
	LongRaw = $('input[name="lon_value"]').val();
	WGSHeightRaw = $('input[name="alt_value"]').val();
// String to Float	
	WGSHeight = parseFloat(WGSHeightRaw);
	Latitude = parseFloat(LatRaw);
	Longitude = parseFloat(LongRaw);

	MajorAxis = 6378137.0 + WGSHeight;
	MinorAxis = 6356752.32 + WGSHeight;
	XEllipse = MajorAxis * dcos(Latitude);
	YEllipse = MinorAxis * dsin(Latitude);
	RadiusEarth = Math.sqrt((XEllipse * XEllipse) + (YEllipse * YEllipse));

// Calculate sineAltitude
// Utilise code from the following website
// Reference: http://www.bogan.ca/astro/telescopes/coodcvtn.html
// Reference: http://www.dept.aoe.vt.edu/~lutze/AOE4134/13LocalSiderealTime.pdf

//	JD = round(days + 2451545.0, 0);
	JD = round(days + 2451545, 0);
	UT = h + min/60;
	JuDay = JD + UT/24;
//	Tut1 = (JD- 2451545.0)/36525;
	Tut1 = JD / 36525;
	UTSidTime = 100.4606184 + 36000.77005361*Tut1 + 0.00038793*Tut1*Tut1;
	UTSidModTime = UTSidTime%360;
	LocSidTime = UTSidModTime + Longitude;
	LocHourAngleMoon = LocSidTime - Ra2;
	LocHourAngleSun = LocSidTime - Ra1;
	SineAltitudeMoon = dsin(Latitude)*dsin(Dec2) + dcos(Latitude)*dcos(Dec2)*dcos(LocHourAngleMoon);
	SineAltitudeSun = dsin(Latitude)*dsin(Dec1) + dcos(Latitude)*dcos(Dec1)*dcos(LocHourAngleSun);
// Calculate Moon Gravity
// 

	EARTHTOMOON = 384400000;
	MASSMOON = 7.3477e+22;
	GCONST = 6.674e-11;
// Use cosine rule to calculate distance to the moon
// EARTHTOMOON^2 = RadiusEarth^2 + DISTANCEMOON^2 - 2*RadiusEarth*DISTANCEMOON*COS(90 + Altitude)
// cos(90+Altitude) = -sin(Altitude)
// Rearrange Equation
// Use Quadratic Formula
// ax^2 + bx + c = 0 

	aQuad = 1;
	bQuad = 2*RadiusEarth*SineAltitudeMoon;
	cQuad = RadiusEarth*RadiusEarth - EARTHTOMOON*EARTHTOMOON;
	DistMoon = (- bQuad + Math.sqrt(bQuad*bQuad - 4*aQuad*cQuad))/(2*aQuad);
// alternative solution is negative and therefore cannot be a solution (i.e. not a distance)
	GrossForceMoon = GCONST * MASSMOON / (DistMoon*DistMoon);

// note that the minus is to match the direction of gravity

// Calculate Sun's Gravity
	EARTHTOSUN = 149600000000;
	MASSSUN = 1.98855e+30;
	GrossForceSun = GCONST * MASSSUN / (EARTHTOSUN*EARTHTOSUN);
	

	// OUTPUTS SUN & MOON GRAVITY
	NormalForceMoon = - GrossForceMoon*SineAltitudeMoon;
	NormalForceSun = - GrossForceSun*SineAltitudeSun;
	

// FORM OUTPUTS SUN & MOON GRAVITY
	form.SunGravity.value = round(NormalForceSun, 6);
	form.MoonGravity.value = round(NormalForceMoon, 6);
   }

//
// this is the usual days since J2000 function
//

   function day2000(y, m, d, h) {
   var d1, b, c, greg;
   greg = y*10000 + m*100 + d;
    if (m == 1 || m == 2) {
	   y = y - 1;
	   m = m + 12;
     }
//  reverts to Julian calendar before 4th Oct 1582
//  no good for UK, America or Sweeden!

   if (greg > 15821004) {
       a = Math.floor(y/100);
       b = 2 - a  + Math.floor(a/4) }
   else {
       b = 0;
     }
   c = Math.floor(365.25 * y);
   d1 = Math.floor(30.6001 * (m + 1));
   return (b + c + d1 -730550.5 + d + h/24);
    } 

//
//	Leap year detecting function (gregorian calendar)
//  returns 1 for leap year and 0 for non-leap year
//

function isleap(y) {
	var a;
//	assume not a leap year...
	a = 0;
//	...flag leap year candidates...
	if (y % 4 == 0) a = 1;
//	...if year is a century year then not leap...
	if (y % 100 ==0 ) a = 0;
//	...except if century year divisible by 400...
	if (y % 400 == 0) a = 1;
//	...and so done according to Gregory's wishes 
	return a;
	}

//
//	Month and day number checking function
//	This will work OK for Julian or Gregorian
//	providing isleap() is defined appropriately
//	Returns 1 if Month and Day combination OK,
//	and 0 if month and day combination impossible
//
function goodmonthday(y, m, d) {
	var a, leap;
	leap = isleap(y); 
//	assume OK
	a = 1;
//	first deal with zero day number!
	if (d == 0) a = 0;
//	Sort Feburary next
	if ((m==2) && (leap ==1) && (d > 29)) a= 0;
	if ((m==2) && (d > 28) && (leap ==0))	a = 0;
//	then the rest of the months - 30 days...
	if(((m==4) || (m == 6) || (m == 9) || (m==11)) && d > 30) a = 0;
//	...31 days...	
	if (d > 31) a = 0;
//	...and so done
	return a;
	}	

//
// Trigonometric functions working in degrees - this just
// makes implementing the formulas in books easier at the
// cost of some wasted multiplications.
// The 'range' function brings angles into range 0 to 360,
// and an atan2(x,y) function returns arctan in correct
// quadrant. ipart(x) returns smallest integer nearest zero
//

function dsin(x) {
	return Math.sin(Math.PI / 180 * x)
	}

function dcos(x) {
	return Math.cos(Math.PI / 180 * x)
	}

function dtan(x) {
	return Math.tan(Math.PI / 180 * x)
	}

function dasin(x) {
	return 180/ Math.PI * Math.asin(x)
	}

function dacos(x) {
	return 180/ Math.PI * Math.acos(x)
	}

function datan(x) {
	return 180/ Math.PI * Math.atan(x)
	}

function datan2(y, x) {
	var a;
	if ((x == 0) && (y == 0)) {
		return 0;
		}
	else	{
		a = datan(y / x);
		if (x < 0) {
			a = a + 180; 
			}
		if (y < 0 && x > 0) {
			a = a + 360;
			}
		return a;
		}
	}

function ipart(x) {
	var a;
	if (x> 0) {
	    a = Math.floor(x);
		}
	else {
		a = Math.ceil(x);
		}	
	return a;
	}

function range(x) {
	var a, b
	b = x / 360;
	a = 360 * (b - ipart(b));
	if (a  < 0 ) {
		a = a + 360
		}
	return a
	}

//
// round rounds the number num to dp decimal places
// the second line is some C like jiggery pokery I
// found in an O'Reilly book which means if dp is null
// you get 2 decimal places.
//
   function round(num, dp) {
//   dp = (!dp ? 2: dp);
   return Math.round (num * Math.pow(10, dp)) / Math.pow(10, dp);
    }


   // -->
		
		
		
		
		
		
    </script>
  </head>

  <body onload="initialize()">
	
    <section>
	    <div id="one">
	    	
	    	<g:form controller="gravity" action="index">
			    <table>
			      <tbody>                                         
			        <tr class="prop">
			          <td valign="top" class="name"><label for="iata">Lat:</label></td>
			          <td valign="top" >
			              <g:select name="lat1" id="lat_direction" from="${['North','South']}" value="${lat_direction }" /><br/>
			              <g:select name="lat2" id="lat_value" from="${0..80}" value="${lat_value }" />	
			              <input name="lat_value" value="${lat }" type="hidden">	              
			          </td>
			        </tr> 
			        <tr class="prop">
			          <td valign="top" class="name"><label for="city">Lon:</label></td>
			          <td valign="top" >
			              <g:select name="lon1" id="lon_direction" from="${['East','West']}" value="${lon_direction }" /><br/>
			              <g:select name="lon2" id="lon_value" from="${0..170}" value="${lon_value }" />
			              <input name="lon_value" value="${lon }" type="hidden">
			          </td>
			        </tr>
			        <tr class="prop">
			          <td valign="top" class="name"><label for="city">Alt:</label></td>
			          <td valign="top" >
			              <input name="alt_value" value="${alt }">
			          </td>
			        </tr>
			        <tr class="prop">
			          <td valign="top" class="name"><label for="city"> </label></td>
			        </tr>         
			        <tr class="prop">
			          <td valign="top" class="name"><label for="city">Gravity:</label></td>
			          <td valign="top" >
			              ${gravity }
			          </td>
			        </tr>                                                 
			      </tbody>
			    </table>
			    <br/>
			    <input name="ne_lat" value="" type="hidden">
			    <input name="ne_lon" value="" type="hidden">
			    <input name="sw_lat" value="" type="hidden">
			    <input name="sw_lon" value="" type="hidden">
			    <input name="nw_lat" value="" type="hidden">
			    <input name="nw_lon" value="" type="hidden">
			    <input name="se_lat" value="" type="hidden">
			    <input name="se_lon" value="" type="hidden">
			    <br/>
			    <input value="Get Gravity" type="submit">
			</g:form>
	    </div>
	    <div id="map-canvas"></div>
	</section>    

  </body>
</html>

