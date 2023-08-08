# ExploreASL Report Config Settings

## Introduction

This is an example of how to use the function `[x, config] =  xASL_qc_GenerateReport(x[, subject])` to generate a completely custom PDF within ExploreASL. If `subject` isn't given, it will assume `subject = x.SUBJECT`.

### Standard Layout.

The PDF generator prints Page by Page, top to bottom in a line by line approach.  For everything it will only print something if it actually exists in your derivatives directory. The standard layout of JSON-objects as follows. On the top level there is one settings meta-object where you can define settings that apply globally to your Report. Other than that there is one other meta-object that define the pages. 

```json
{   
    "Global_Settings" : {
        "type" : "settings",
    }, 
    "Page_one" : {
        "type" : "page", 
        "pageIdentifier" : "1", 
        "content" : {[
            {   "type" : "text",
                "text" : "something"},
            {   "type" : "image",
                "path" : "<path to image>"}
        ]}
    }, 
    "Page_two" : {
        "type" : "page", 
        "pageIdentifier" : "2", 
        "content" : {[
            {   "type" : "text",
                "text" : "something"},
            {   "type" : "image",
                "path" : "<path to image>"}
        ]}
    }
}
```
For all objects in the content array a "type" is needed that defines what kind of object it pertains. Currently the following objects are recognized : `settings, page, image, scan, text, block`. The options for the contents of each object will be discussed below.

## Meta Objects
### Page
The `page` object is Meta Object used to create a new page of whatever you want to print. Without page elements nothing will be printed. A page element cannot contain other pages itself.

```json
{
    "type" : "page", 
    "pageIdentifier" : "numerical - REQUIRED", 
    "content" : {[
            {   "type" : "text",
                "text" : "something"},
            {   "type" : "image",
                "path" : "<path to image>"}
    ]}
}
```

`type` is used to define the object as a page object. 

`pageIdentifier` defines the page name of the final pdf. 

`content` is an array of Content objects defined below, printed in order top to bottom. 

### Settings
The `settings` object can be called as a Meta Object, to set print settings globally. It can also be called as a Content object to change the settings of a page or block. It can also be called withing a Content object to change the settings of just that specific object.  
This object has the following properties and value types : 

```json
{
    "type" : "settings", 
    "fontSize" : "numercal", 
    "fontName" : "string", 
    "fontWeight" : "string", 
    "color" : "string", 
    "axesVisible" : "string", 
    "HorizontalAlignment" : "string", 
    "lineSpacing" : "numerical"
}
```

`type` is used to define the object as a settings object.

 `fontSize, lineSpacing` need numerical values and determine fontSize in px and lineSpacing in % of the full page.

`fontName` can be any MATLAB accepted font. 

`fontWeight` can be `normal, bold` and determines bolding of text. 

`color` can take any matlab supported color format:
- long names `red, blue, green ...` 
- short names `r, g, b, k, w ...` 
- hex values `#FF000`

`HorizontalAlignment` can be set to `left, right, center`. 

`axesVisible` can be use to see where your text will be printed for debugging purposes, accepted values are `off, on`. 

If a property is not defined, the (previously defined) default will be used.



## Content Objects
### Text
```json
 {
    "type" : "text", 
    "settings" : {
        "fontSize" : "numercal", 
        "fontWeight" : "string"
    }, 
    "text" : "string REQUIRED"
}
```

`type` is used to define the object as a printable text object.  

`settings` can be used to change the print settings for this line only. 

`text` contains the string to be printed.

### Image
```json
{
    "type" : "image", , 
    "absolutePath" : "string Semi-REQUIRED", 
    "popPath" : "string Semi-REQUIRED",
    "subjPath" : "string Semi-REQUIRED",
    "position" : "[x y] REQUIRED", 
    "size" : "[x y] REQUIRED"    
}
```
`type` is used to define the object as a printable image object. 

`absolutePath`, `popPath` and `subjPath` define the location of the image you want to print to page. 

**Exactly _one_ of these three must be given.**
- Use `absolutePath` if the file is on a fixed location, e.g. `~/images/example.jpg`). 
- Use `popPath` if the file is located in the population section e.g. `M0/M0check.jpg`).
- Use `subjPath` if the file is located in the population section e.g. `M0/M0check.jpg`).

For all three of these you can use wildcards (formatted as`<string>`)to replace strings from the x-struct (eg. `M0Check/Cor_M0_<SUBJECT>_<SESSION>` would translate to `M0Check/Cor_M0_1_ASL_1` if `x.SUBJECT = '1'` & `x.SESSION ='ASL_1'`.

`position` needs a 2 element array with values between 0 and 1 of print location relative to the page. (eg. `[0.5 0.5]` would be the center of the page). 

`size` needs a similar 2 element array defining the size of the printed image.

### Scan
```json
{
    "type" : "scan", 
    "name" : "string REQUIRED", 
    "position" : "[x y] REQUIRED", 
    "size" : "[x y] REQUIRED", 
    "slice" : {
        "TraSlice" : "[n1 n2 n3 ...] REQUIRED", 
        "CorSlice" : "[n1 n2 n3 ...] REQUIRED", 
        "SagSlice" : "[n1 n2 n3 ...] REQUIRED"}, 
    "overlay" : {
        "string" : {
            "contour" : "bool", 
            "opacity" : "numerical", 
            "color" : "string"
        }
    }
}
```
`type` is used to define the object as a printable scan object. 

`name` is used to identify the scan to be used, any NIfTI file in the x.P module can be used. (eg. `Pop_Path_rT1` for the rT1 scan in the population module)

`slice` has three elements, one for each scan orientation, you can define any amount of slices to be printed, `"TraSlice" : "[25 75], CorSlice" : "[50], "SagSlice" : "[]"` would print the transversal slices at 25% and 75% position, the central coronal slice and no saggital slices for exmple.

`position` needs a 2 element array with values between 0 and 1 of print location relative to the page. (eg. `[0.5 0.5]` would be the center of the page). `size` needs a similar 2 element array defining the size of the printed image.

`overlay` is an optional property to be used to overlay any other scan over the original. The `string` used as property can be any NIfTI file in the population module, just like the `name` property above. `contour` can take values `true, false`. `opacity` takes a value between 0 and 1. `color` uses any MATLAB defined color just like the settings object defined above. Example overlay : `"SliceGradient" : { "contour" : "true", "opacity" : "0.3", "color" : "r" }` 

### Block
```json
{
    "type" : "block", 
    "position" : "[0.2 0.11]", 
    "size" : "[0.8 0.33]", 
    "settings": {
	      "type" : "settings",
			  "<any setting seen in the previous section>" : "parameter"
},
"content" : {[
            {   "type" : "text",
                "text" : "something"},
            {   "type" : "image",
                "path" : "<path to image>"}
   ]}
}
```

A block object creates a subfigure, saves it as an image, and prints the image in your predefined position, you can print anything you want here. works identical to a `page` object in that you can add any number of printable objects or settings objects. 

`type` is used to define the object as a printable block object. 

`position` needs a 2 element array with values between 0 and 1 of print location relative to the page. (eg. `[0.5 0.5]` would be the center of the page). `size` needs a similar 2 element array defining the size of the printed image.

`content` is an array of content objects, printed in order top to bottom. Functions identical to the implementation of `content` in `page`. As `block` is a content object itself you can recursively put blocks inside one another.
