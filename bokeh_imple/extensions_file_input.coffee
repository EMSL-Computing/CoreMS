import * as p from "core/properties"
import {WidgetBox, WidgetBoxView} from "models/layouts/widget_box"

export class FileInputView extends WidgetBoxView

  initialize: (options) ->
    super(options)
    input = document.createElement("input")
    input.type = "file"
    input.accept = @model.accept
    input.id = @model.id
    input.style = "width:" + @model.width + "px"
    input.onchange = () =>
      @model.value = input.value
      @model.file_name = input.files[0].name
      @file_handler(input)
    @el.appendChild(input)

  file_handler: (input) ->
    file = input.files[0]
    opts =  
      header: true,
      dynamicTyping: true,
      delimiter: ",",
      newline: "\r\n",
      complete: (results) =>
        input.data = results.data
        @.model.data = results.data
    Papa.parse(file, opts)


export class FileInput extends WidgetBox
  default_view: FileInputView
  type: "FileInput"
  @define {
    value: [ p.String ]
    file_name: [ p.String ]
    accept: [ p.String ]
    data : [ p.Array ]
  }