<OpenGeoSysGLI>
    <name>cuboid_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="-60" />
        <point id="1" x="0" y="0" z="0" />
        <point id="2" x="0" y="2" z="0" />
        <point id="3" x="0" y="2" z="-60" />
		<point id="4" x="2" y="0" z="-60" />
		<point id="5" x="2" y="0" z="0" />
		<point id="6" x="2" y="2" z="0" />
		<point id="7" x="2" y="2" z="-60" />
    </points>

    <surfaces>
        <surface id="0" name="left">
            <element p1="0" p2="1" p3="2" />
            <element p1="0" p2="3" p3="2" />
        </surface>
        <surface id="1" name="right">
            <element p1="4" p2="6" p3="5" />
            <element p1="4" p2="6" p3="7" />
        </surface>
        <surface id="2" name="top">
            <element p1="1" p2="2" p3="5" />
            <element p1="5" p2="2" p3="6" />
        </surface>
        <surface id="3" name="bottom">
            <element p1="0" p2="3" p3="4" />
            <element p1="4" p2="3" p3="7" />
        </surface>
        <surface id="4" name="front">
            <element p1="0" p2="1" p3="4" />
            <element p1="4" p2="1" p3="5" />
        </surface>
        <surface id="5" name="back">
            <element p1="2" p2="3" p3="6" />
            <element p1="6" p2="3" p3="7" />
        </surface>
    </surfaces>
</OpenGeoSysGLI>