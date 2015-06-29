var greatFix = function(){
  var data = $$('form').retrieve('gold')[0].options.unitData;
  $$('.checkboxes').each(function(cb){
    cb.getElements('input').each(function(input){
      if (data[input.get('class').split(" ")[0]] == input.value) {
        input.checked = true;
        input.fireEvent('change')
      }
    })
  })
}

if(_cf_cml.digging_gold){
  greatFix.delay(100)
}
