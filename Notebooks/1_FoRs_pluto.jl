### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ e6368433-ce9b-4a90-b226-32ee77ac3e37
using Geodesy

# ╔═╡ 501a8fa3-abf8-4858-8723-5d8be0ac4f78
using CoordinateTransformations, Rotations, StaticArrays

# ╔═╡ 2de58f20-7bab-11ec-1dde-dbf8aef620e9
md"""
We illustrate here some manipulations of the standard frames of reference (FoRs) and coordinate computations in these frames. We illustrate in particular some of the functionality of the [Geodesy.jl](https://github.com/JuliaGeo/Geodesy.jl) package.
"""

# ╔═╡ cb7b0b5d-b901-44b4-88a5-e3314f16ec73
md"""
# ECEF FoR and LLA coordinates
"""

# ╔═╡ 94dd5109-3425-4c88-9a1a-fdab06320e6c
md"""
The ECEF frame comes with a natural rectangular coordinate system. But positions on the surface of the Earth are more naturally specified in the ellipsoïdal coordinate system: latiture, longitude, altitude (LLA).
"""

# ╔═╡ 9a5f85b4-579b-439b-a4b0-1dbf7b8160d8
@doc ECEF

# ╔═╡ b1972673-2247-4ea6-af9c-9abcbb9c4f8a
@doc LLA

# ╔═╡ 4439c783-8947-4efd-b61a-cf0c46b6f6af
poly_lla = LLA(45.50439, -73.61288, 159.0)

# ╔═╡ 9f423556-cae1-41b0-be2f-f1bf966c9770
md"""
Conversion from LLA to ECEF rectangular coordinates requires the specification of a **reference ellipsoid** provided by a datum, such as WGS 84. Recall that WGS 84 is the datum used by the Global Positioning System.
"""

# ╔═╡ 7a6b7a99-fc23-458a-85b7-104a3169c20f
poly_ecef = ECEF(poly_lla, wgs84)

# ╔═╡ 64f119df-b5c1-4512-9b82-9eb459c9025b
md"""
Here is the (slightly different) result for the [OGSB 36](https://en.wikipedia.org/wiki/Ordnance_Survey_National_Grid) datum.
"""

# ╔═╡ c3d1cb89-c24a-4fc6-9963-ff8c0be81a64
poly_ecef2 = ECEF(poly_lla, osgb36)

# ╔═╡ 35d466ec-3f16-4404-8e5f-d6ca8158fb88
md"""
We can define a coordinate transformation as follows:
"""

# ╔═╡ ecb9a286-a8ed-4aa9-b139-85f24bf9d7c1
T₁ = ECEFfromLLA(wgs84)

# ╔═╡ 7ac54799-bad1-43b4-978c-9dd01e8169a0
# Example
T₁(poly_lla)

# ╔═╡ 2c829b32-2d63-4021-ac75-b00c4eb6555c
md"""
The [Universal Transverse Mercator](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system) is another coordinate system used to identify locations on the surface of the Earth. See also the [Universal Polar Stereographic](https://en.wikipedia.org/wiki/Universal_polar_stereographic_coordinate_system) coordinate system.
"""

# ╔═╡ 52675ea4-8450-4b69-b9c5-36977700bb2b
poly_utmz = UTMZ(poly_lla, wgs84)

# ╔═╡ fa07ff7b-476d-44b5-8947-0b2009787f0f
md"""
# Local Tangent Frames
"""

# ╔═╡ 784eed52-d0d0-4c0a-b3c7-c1d1622ce8c9
md"""
These frames are NED (North-East-Down) or ENU (East-North-Up) frames. The Geodesy.jl package supports ENU.
"""

# ╔═╡ 070cb0d9-9a30-4df6-8703-16f6a9639d37
@doc ENU

# ╔═╡ a5a4763e-0ca0-4348-ab90-acacd76d7b67
md"""
Coordinates of McGill in ENU frame centered at Polytechnique:
"""

# ╔═╡ a2d3986f-d214-4c46-b7cc-d10fd7818f1d
begin
mcgill_lla = LLA(45.5047847,-73.5771511,47.9)
mcgill_polyENU = ENU(mcgill_lla, poly_lla, wgs84)  
end

# ╔═╡ 99b7eabf-f1ef-4ae6-99a0-b018f76710ad
begin
using LinearAlgebra
norm(mcgill_polyENU)
end

# ╔═╡ d8c37bd5-dba7-4290-9800-ba727c0bb632
md"""
Distance between McGill and Polytechnique:
"""

# ╔═╡ bc4f6aec-ef8f-44fe-b169-a89327b72928
md"""
What is the z-component of McGill in Polytechnique's ENU frame that is due solely to the curvature of the Earth?
"""

# ╔═╡ 9d0b120c-1408-4842-8cd6-7dab106207cc
abs(mcgill_polyENU.u) - abs(poly_lla.alt - mcgill_lla.alt)

# ╔═╡ b318a3ee-d329-4947-af77-ee2dd9ed4504
md"""
We can cefine a coordinate transformation from LLA to the ENU frame centered at Polytechnique:
"""

# ╔═╡ a040572d-e897-4ca8-9443-9bf1ce0fe6d3
T₂ = ENUfromLLA(poly_lla, wgs84)

# ╔═╡ ddb23512-97d0-4e59-b50b-f9ecc27c8706
# Example: Quebec city in Poly ENU
begin
	quebec_lla = LLA(46.829853, -71.254028, 74)
	quebec_polyENU = T₂(quebec_lla)
end

# ╔═╡ 2c75e406-97b7-40ff-a02f-ffffdf2349ad
md"""
Let's check the earlier result about the influence of the Earth curvature, but now using the functions of the toolbox. First, we define an ENU frame centered at Polytechnique but at altitude 0.
"""

# ╔═╡ 2529c01c-482a-4695-b2a7-56ec7dbfc4dc
T₃ = ENUfromLLA(LLA(poly_lla.lat, poly_lla.lon, 0.0), wgs84) 

# ╔═╡ 28f7d9e0-1912-4229-bd65-08702c221b9b
md"""
Here is how much below the z-plane of poly_ENU is McGill:
"""

# ╔═╡ 65c57035-868d-4307-a967-7a140f4427a4
T₃(LLA(mcgill_lla.lat, mcgill_lla.lon, 0.0)).u   

# ╔═╡ 22b66270-e441-4132-8927-15a613b7db28
md"""
# Instrument Frames and Coordinate Transformations
"""

# ╔═╡ 233f9295-1c39-44d9-ae7e-94aa2789422d
md"""
The frames of reference above and coordinate transformations are compatible with the [CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl) package. Hence, they define coordinate transformations that can be composed with custom ones.
"""

# ╔═╡ 0b294bd7-d490-493b-a6c9-c18b2962e499
md"""
Let us assume for instance that we have an aerial vehicle with its own FoR (forward, right, down). On it is mounted a sensor that can measure the spherical coordinates (azimuth, elevation, distance) of an object in this FoR (for example, a LiDAR). We want to define a transformation  pipeline that takes a sensor measurement and the knowledge of the vehicle LLA coordinates and vehicle orientation, to directly return the *ECEF rectangular coordinates* of the measurement.
"""

# ╔═╡ 203fb0d6-bd15-49bc-bc2d-2ca673c4c57e
md"""
We have a sensor frame with a certain fixed orientation and position in the vehicle frame:
"""

# ╔═╡ 1a31b1b7-cb4b-4340-ba4f-95cccb5f3097
begin
	# Yaw, pitch, roll angles for the sensor
	ϕsensor, θsensor, ψsensor = π/3, π/5, -π/4   
	# Coordinates of the center of the sensor frame in the vehicle frame
	pSensor = SVector{3}([1,2,0.5])	
end

# ╔═╡ 94425731-4016-42d7-96b4-1cb7f3da8329
md"""
We define a coordinate transformation from the sensor to the navigation frame. This is also the rigid transformation transporting the navigation frame to the sensor frame.
"""

# ╔═╡ 421b9657-93db-47fb-a5ce-052d4a4f5d9a
begin
	rSensor = RotZYX(ϕsensor, θsensor, ψsensor)
	sensorToVehicle = AffineMap(rSensor,pSensor)
end

# ╔═╡ f09df84e-934f-4948-a21d-ba6dd5a74437
md"""
So, here are the coordinates in the vehicle frame of a point 10m on the z axis of the sensor frame.
"""

# ╔═╡ e9e816af-4c5a-452e-97f9-a19b5ae9780d
sensorToVehicle([0,0,10])

# ╔═╡ 4fd92dbb-b768-47e9-a9c0-930053f4508e
md"""
The sensor takes measurements in spherical coordinates. Consider the spherical coordinates defined in CoordinateTransformations.jl. In general, we would have to make sure our sensor conventions are the same.
"""

# ╔═╡ 349569d0-4471-4845-b59b-fba1e207593d
md"""
Spherical coordinates from CoordinateTransformations.jl:
"""

# ╔═╡ b316eb2b-3f9e-4fc0-98a9-786481d5cf4a
@doc Spherical

# ╔═╡ be6a9a6b-e8ad-4817-9111-9552f1facb43
targetMeasurement = Spherical(1,π/3,π/4)

# ╔═╡ a18b2b68-d8d5-48d9-b91c-08c921cb79ef
CartesianFromSpherical()(targetMeasurement)

# ╔═╡ fc9cef22-4aa6-4f51-8cfd-1dbe631616d1
md"""
The transformation we are looking for works as follows. We first transform the measurements from spherical coordinate to rectangular coordinates, in the sensor frame. We then change coordinates to the vehicle frame.
"""

# ╔═╡ a4a91ab5-455e-4c9e-aa48-7e66bc384890
measurementToVehicle = sensorToVehicle ∘ CartesianFromSpherical()

# ╔═╡ 6d863bf0-9a7c-4c16-b194-629cb4168190
measurementToVehicle(targetMeasurement)

# ╔═╡ a78ec4dd-eefe-4e81-8caf-52854c914849
md"""
At this point, we want to change coordinates from the vehicle frame to the geographic frame, an ENU frame centered at the projection of the vehicle frame on the ellipsoid. This requires the orientation of the vehicle as well as its altitude.
"""

# ╔═╡ 22de5670-607c-496f-80f8-28483fd33bf9
begin
	yawVehi, pitchVehi, rollVehi = π/6, π/20, π/10
	rVehicle = RotZYX(yawVehi, pitchVehi, rollVehi)
end

# ╔═╡ cb766458-5a31-4826-b599-afb10e1ccd39
pVehicle_lla = LLA(45.50439, -73.61288, 159.0)  # at Polytechnique

# ╔═╡ bcb6e38e-7ab5-4a79-b822-2622b887e131
vehiToGeographic = AffineMap(RotZYX(yawVehi, pitchVehi, rollVehi),SVector{3}([0,0,pVehicle_lla.alt]))

# ╔═╡ 34875fd5-da63-4060-96ca-9fbcf065eb11
measurement_enu = (vehiToGeographic ∘ measurementToVehicle)(targetMeasurement)

# ╔═╡ e286174f-57f0-48c9-b92b-cef2455ed6cf
md"""
Finally, to obtain the coordinates of the measurement in the rectangular ECEF coordinates, this requires the LLA coordinates of the center of the last ENU frame (position of the vehicle).
"""

# ╔═╡ eedf6d96-1e76-46a9-85b3-688b7663f4e0
ECEFfromENU(pVehicle_lla,wgs84)(ENU(measurement_enu))

# ╔═╡ 100911ab-b58e-42ac-9624-dce521424eea
function sensorMeasToECEF(r, θ, ϕ, sensorOrientation, sensorPosition, vehicle_lla, vehicle_orientation; datum=wgs84)
	sensorToVehicle = AffineMap(sensorOrientation, sensorPosition)

	# ENU frame here is at altitude 0
	vehicleToGeographic = AffineMap(vehicle_orientation, SVector{3}([0,0, vehicle_lla.alt]))
	
	measurementToGeographic = vehicleToGeographic ∘ sensorToVehicle ∘ CartesianFromSpherical()

	measurement_enu = ENU(measurementToGeographic(Spherical(r, θ, ϕ)))
	return ECEFfromENU(vehicle_lla, datum)(ENU(measurement_enu))
end

# ╔═╡ 2b0f31b2-9f20-495a-9588-dbb692105100
md"""
Enter the sensor measurements:
"""

# ╔═╡ 894e71a2-b73b-4a13-85c0-d1c54b64f21d
r, ϕ, ψ = 100, π/3, π/4

# ╔═╡ ff8d7258-3d61-490f-90f7-ad44884940d5
md"""
ECEF coordinates of the sensor measurements: 

$(sensorMeasToECEF(r, ϕ, ψ, rSensor, pSensor, pVehicle_lla, rVehicle))
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
Geodesy = "0ef565a4-170c-5f04-8de2-149903a85f3d"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
CoordinateTransformations = "~0.6.2"
Geodesy = "~1.0.1"
Rotations = "~1.1.1"
StaticArrays = "~1.3.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "54fc4400de6e5c3e27be6047da2ef6ba355511f8"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.Geodesy]]
deps = ["CoordinateTransformations", "Dates", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "2341a0d40d1f96db72275d7ae97ad5600778f137"
uuid = "0ef565a4-170c-5f04-8de2-149903a85f3d"
version = "1.0.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra"]
git-tree-sha1 = "adf644ef95a5e26c8774890a509a55b7791a139f"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "2fa87d198bc5356c649b92109ed3ce46ee1eb89d"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.1.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2884859916598f974858ff01df7dfc6c708dd895"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─2de58f20-7bab-11ec-1dde-dbf8aef620e9
# ╠═e6368433-ce9b-4a90-b226-32ee77ac3e37
# ╟─cb7b0b5d-b901-44b4-88a5-e3314f16ec73
# ╟─94dd5109-3425-4c88-9a1a-fdab06320e6c
# ╟─9a5f85b4-579b-439b-a4b0-1dbf7b8160d8
# ╟─b1972673-2247-4ea6-af9c-9abcbb9c4f8a
# ╠═4439c783-8947-4efd-b61a-cf0c46b6f6af
# ╟─9f423556-cae1-41b0-be2f-f1bf966c9770
# ╠═7a6b7a99-fc23-458a-85b7-104a3169c20f
# ╟─64f119df-b5c1-4512-9b82-9eb459c9025b
# ╠═c3d1cb89-c24a-4fc6-9963-ff8c0be81a64
# ╟─35d466ec-3f16-4404-8e5f-d6ca8158fb88
# ╠═ecb9a286-a8ed-4aa9-b139-85f24bf9d7c1
# ╠═7ac54799-bad1-43b4-978c-9dd01e8169a0
# ╟─2c829b32-2d63-4021-ac75-b00c4eb6555c
# ╠═52675ea4-8450-4b69-b9c5-36977700bb2b
# ╟─fa07ff7b-476d-44b5-8947-0b2009787f0f
# ╟─784eed52-d0d0-4c0a-b3c7-c1d1622ce8c9
# ╟─070cb0d9-9a30-4df6-8703-16f6a9639d37
# ╟─a5a4763e-0ca0-4348-ab90-acacd76d7b67
# ╠═a2d3986f-d214-4c46-b7cc-d10fd7818f1d
# ╟─d8c37bd5-dba7-4290-9800-ba727c0bb632
# ╠═99b7eabf-f1ef-4ae6-99a0-b018f76710ad
# ╟─bc4f6aec-ef8f-44fe-b169-a89327b72928
# ╠═9d0b120c-1408-4842-8cd6-7dab106207cc
# ╟─b318a3ee-d329-4947-af77-ee2dd9ed4504
# ╠═a040572d-e897-4ca8-9443-9bf1ce0fe6d3
# ╠═ddb23512-97d0-4e59-b50b-f9ecc27c8706
# ╟─2c75e406-97b7-40ff-a02f-ffffdf2349ad
# ╠═2529c01c-482a-4695-b2a7-56ec7dbfc4dc
# ╟─28f7d9e0-1912-4229-bd65-08702c221b9b
# ╠═65c57035-868d-4307-a967-7a140f4427a4
# ╟─22b66270-e441-4132-8927-15a613b7db28
# ╟─233f9295-1c39-44d9-ae7e-94aa2789422d
# ╟─0b294bd7-d490-493b-a6c9-c18b2962e499
# ╠═501a8fa3-abf8-4858-8723-5d8be0ac4f78
# ╟─203fb0d6-bd15-49bc-bc2d-2ca673c4c57e
# ╠═1a31b1b7-cb4b-4340-ba4f-95cccb5f3097
# ╟─94425731-4016-42d7-96b4-1cb7f3da8329
# ╠═421b9657-93db-47fb-a5ce-052d4a4f5d9a
# ╟─f09df84e-934f-4948-a21d-ba6dd5a74437
# ╠═e9e816af-4c5a-452e-97f9-a19b5ae9780d
# ╟─4fd92dbb-b768-47e9-a9c0-930053f4508e
# ╟─349569d0-4471-4845-b59b-fba1e207593d
# ╟─b316eb2b-3f9e-4fc0-98a9-786481d5cf4a
# ╠═be6a9a6b-e8ad-4817-9111-9552f1facb43
# ╠═a18b2b68-d8d5-48d9-b91c-08c921cb79ef
# ╟─fc9cef22-4aa6-4f51-8cfd-1dbe631616d1
# ╠═a4a91ab5-455e-4c9e-aa48-7e66bc384890
# ╠═6d863bf0-9a7c-4c16-b194-629cb4168190
# ╟─a78ec4dd-eefe-4e81-8caf-52854c914849
# ╠═22de5670-607c-496f-80f8-28483fd33bf9
# ╠═cb766458-5a31-4826-b599-afb10e1ccd39
# ╠═bcb6e38e-7ab5-4a79-b822-2622b887e131
# ╠═34875fd5-da63-4060-96ca-9fbcf065eb11
# ╟─e286174f-57f0-48c9-b92b-cef2455ed6cf
# ╠═eedf6d96-1e76-46a9-85b3-688b7663f4e0
# ╠═100911ab-b58e-42ac-9624-dce521424eea
# ╟─2b0f31b2-9f20-495a-9588-dbb692105100
# ╠═894e71a2-b73b-4a13-85c0-d1c54b64f21d
# ╟─ff8d7258-3d61-490f-90f7-ad44884940d5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
