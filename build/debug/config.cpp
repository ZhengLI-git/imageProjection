{
    "xray":{
		"d_f3d": 50,//Distance of ray source (focal point) to 3D object
		"d_3d2d":100,//Distance of 3D object to DRR
		"npixel_x":1024,// number of pixels along X of the 2D DRR image
		"npixel_y":1024,// number of pixels along Y of the 2D DRR image
		"res_x":0.5,// pixel spacing along X of the 2D DRR image [mm]
		"res_y":0.5// pixel spacing along Y of the 2D DRR image [mm]
    },
    "transformation":{
		"tx_min":0,//Translation parameter of the object
		"tx_max":0, 
		"tx_delta":0, 
		"ty_min":0,
		"ty_max":0, 
		"ty_delta":0, 
		"tz_min":0, 
		"tz_max":0, 
		"tz_delta":0,
		
		"rx_min":0, //Rotation around x,y,z axis in degrees 
		"rx_max":0, 
		"rx_delta":0, 
		"ry_min":-30, 
		"ry_max":30, 
		"ry_delta":5, 
		"rz_min":-30, 
		"rz_max":30, 
		"rz_delta":5
	}
}
