It SHOULD build if you navigate to this folder in your command line and run `cmake -S ./ -B ./Build`

In case .obj file imports for the horse and statue are not linked properly when you open Houdini, the .obj files have been included in the submission.

Included in the submission is /dso/ShatterFX.dll

If you have Houdini version 20.5.487, this .dll should be able to be added to your Houdini /dso/ folder, enabling the use of our custom nodes in Houdini. If this doesn't work or you cannot access this version of Houdini, you may need to recompile our solution with your local version.