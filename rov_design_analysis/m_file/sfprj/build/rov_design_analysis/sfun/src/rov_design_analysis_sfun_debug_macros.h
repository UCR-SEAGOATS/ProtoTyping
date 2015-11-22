#ifndef __SF_DEBUG_MACROS_H__
#define __SF_DEBUG_MACROS_H__

#define _SFD_MACHINE_CALL(v1,v2,v3) sf_debug_call(_rov_design_analysisMachineNumber_,UNREASONABLE_NUMBER,UNREASONABLE_NUMBER,MACHINE_OBJECT,v1,v2,v3,(unsigned int) _sfEvent_,-1,NULL,_sfTime_,1)
#define _SFD_ME_CALL(v2,v3) _SFD_MACHINE_CALL(EVENT_OBJECT,v2,v3)
#define _SFD_MD_CALL(v2,v3) _SFD_MACHINE_CALL(EVENT_OBJECT,v2,v3)
extern unsigned int _rov_design_analysisMachineNumber_;
#define _SFD_SET_DATA_VALUE_PTR(v1,v2)\
	sf_debug_set_instance_data_value_ptr(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,CHARTINSTANCE_INSTANCENUMBER,v1,(void *)(v2));
#define _SFD_UNSET_DATA_VALUE_PTR(v1)\
	sf_debug_unset_instance_data_value_ptr(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,CHARTINSTANCE_INSTANCENUMBER,v1);
#define _SFD_DATA_RANGE_CHECK_MIN_MAX(dVal,dNum,dMin,dMax)\
                      sf_debug_data_range_error_wrapper_min_max(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             dNum,(double)(dVal),(double)dMin,(double)dMax)
#define _SFD_DATA_RANGE_CHECK_MIN(dVal,dNum,dMin)\
                      sf_debug_data_range_error_wrapper_min(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             dNum,(double)(dVal),(double)dMin)
#define _SFD_DATA_RANGE_CHECK_MAX(dVal,dNum,dMax)\
                      sf_debug_data_range_error_wrapper_max(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             dNum,(double)(dVal),(double)dMax)
#define _SFD_DATA_RANGE_CHECK(dVal,dNum)\
                      sf_debug_data_range_wrapper(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             dNum,(double)(dVal))
#define _SFD_ARRAY_BOUNDS_CHECK(v1,v2,v3,v4,v5) \
                      sf_debug_data_array_bounds_error_check(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),(int)(v2),(int)(v3),(int)(v4),(int)(v5))
#define _SFD_EML_ARRAY_BOUNDS_CHECK(v1,v2,v3,v4,v5) \
                      sf_debug_eml_data_array_bounds_error_check(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),(int)(v2),(int)(v3),(int)(v4),(int)(v5))
#define _SFD_INTEGER_CHECK(v1,v2) \
                      sf_debug_integer_check(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),(double)(v2))
#define _SFD_CAST_TO_UINT8(v1) \
                      sf_debug_cast_to_uint8_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_UINT16(v1) \
                      sf_debug_cast_to_uint16_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_UINT32(v1) \
                      sf_debug_cast_to_uint32_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_INT8(v1) \
                      sf_debug_cast_to_int8_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_INT16(v1) \
                      sf_debug_cast_to_int16_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_INT32(v1) \
                      sf_debug_cast_to_int32_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_CAST_TO_SINGLE(v1) \
                      sf_debug_cast_to_real32_T(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
                                             (v1),0,0)
#define _SFD_TRANSITION_CONFLICT(v1,v2) sf_debug_transition_conflict_error(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
v1,v2)
#define _SFD_CHART_CALL(v1,v2,v3) sf_debug_call(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
CHART_OBJECT,v1,v2,v3,(unsigned int)_sfEvent_,\
0,NULL,_sfTime_,1)
#define _SFD_CC_CALL(v2,v3) _SFD_CHART_CALL(CHART_OBJECT,v2,v3)
#define _SFD_CS_CALL(v2,v3) _SFD_CHART_CALL(STATE_OBJECT,v2,v3)
#define _SFD_CT_CALL(v2,v3) _SFD_CHART_CALL(TRANSITION_OBJECT,v2,v3)
#define _SFD_CE_CALL(v2,v3) _SFD_CHART_CALL(EVENT_OBJECT,v2,v3)
#define _SFD_CD_CALL(v2,v3) _SFD_CHART_CALL(EVENT_OBJECT,v2,v3)
#define _SFD_EML_CALL(v1,v2,v3) sf_debug_call(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
CHART_OBJECT,STATE_OBJECT,v1,v2,(unsigned int)_sfEvent_,\
v3,NULL,_sfTime_,1)
#define _SFD_CHART_COVERAGE_CALL(v1,v2,v3,v4) sf_debug_call(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
CHART_OBJECT,v1,v2,v3,(unsigned int) _sfEvent_,\
v4,NULL,_sfTime_,1)
#define _SFD_CCS_CALL(v2,v3,v4) _SFD_CHART_COVERAGE_CALL(STATE_OBJECT,v2,v3,v4)
#define _SFD_CCT_CALL(v2,v3,v4) _SFD_CHART_COVERAGE_CALL(TRANSITION_OBJECT,v2,v3,v4)
#define _SFD_CCP_CALL(v3,v4,v5) sf_debug_call(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
CHART_OBJECT,TRANSITION_OBJECT,TRANSITION_GUARD_COVERAGE_TAG,v3,(unsigned int) _sfEvent_,\
v4,NULL,_sfTime_,(unsigned int)(v5))
#define _SFD_STATE_TEMPORAL_THRESHOLD(v1,v2,v4) sf_debug_temporal_threshold(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
(unsigned int)(v1),(v2),STATE_OBJECT,(v4))
#define _SFD_TRANS_TEMPORAL_THRESHOLD(v1,v2,v4) sf_debug_temporal_threshold(_rov_design_analysisMachineNumber_,\
CHARTINSTANCE_CHARTNUMBER,\
CHARTINSTANCE_INSTANCENUMBER,\
(unsigned int)(v1),(v2),TRANSITION_OBJECT,(v4))
#define CV_EVAL(v1,v2,v3,v4) cv_eval_point(_rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(boolean_T)(v4))
#define CV_CHART_EVAL(v2,v3,v4) CV_EVAL(CHART_OBJECT,(v2),(v3),(v4))
#define CV_STATE_EVAL(v2,v3,v4) CV_EVAL(STATE_OBJECT,(v2),(v3),(v4))
#define CV_TRANSITION_EVAL(v1,v2) cv_eval_point(_rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  TRANSITION_OBJECT,(v1),0,((v2)!=0))

/* Coverage EML Macros */
#define CV_EML_EVAL(v1,v2,v3,v4) cv_eml_eval(_rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(int)(v4))
#define CV_EML_FCN(v2,v3) CV_EML_EVAL(CV_EML_FCN_CHECK,(v2),(v3),0)
#define CV_EML_IF(v2,v3,v4) CV_EML_EVAL(CV_EML_IF_CHECK,(v2),(v3),(v4))
#define CV_EML_FOR(v2,v3,v4) CV_EML_EVAL(CV_EML_FOR_CHECK,(v2),(v3),(v4))
#define CV_EML_WHILE(v2,v3,v4) CV_EML_EVAL(CV_EML_WHILE_CHECK,(v2),(v3),(v4))
#define CV_EML_SWITCH(v2,v3,v4) CV_EML_EVAL(CV_EML_SWITCH_CHECK,(v2),(v3),(v4))
#define CV_EML_COND(v2,v3,v4) CV_EML_EVAL(CV_EML_COND_CHECK,(v2),(v3),(v4))
#define CV_EML_MCDC(v2,v3,v4) CV_EML_EVAL(CV_EML_MCDC_CHECK,(v2),(v3),(v4))

#define _SFD_CV_INIT_EML(v1,v2,v3,v4,v5,v6,v7,v8) cv_eml_init_script(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5),(v6),(v7),(v8))

#define _SFD_CV_INIT_EML_FCN(v1,v2,v3,v4,v5,v6) cv_eml_init_fcn(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5),(v6))

#define _SFD_CV_INIT_EML_IF(v1,v2,v3,v4,v5,v6) cv_eml_init_if(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5),(v6))

#define _SFD_CV_INIT_EML_FOR(v1,v2,v3,v4,v5) cv_eml_init_for(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5))

#define _SFD_CV_INIT_EML_WHILE(v1,v2,v3,v4,v5) cv_eml_init_while(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5))

#define _SFD_CV_INIT_EML_MCDC(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) cv_eml_init_mcdc(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5),(v6),(v7),(v8),(v9),(v10))

#define _SFD_CV_INIT_EML_SWITCH(v1,v2,v3,v4,v5,v6,v7,v8) cv_eml_init_switch(\
       _rov_design_analysisMachineNumber_,\
		  CHARTINSTANCE_CHARTNUMBER,\
		  CHARTINSTANCE_INSTANCENUMBER,\
		  (v1),(v2),(v3),(v4),(v5),(v6),(v7),(v8))


#define _SFD_SET_DATA_PROPS(dataNumber,dataScope,isInputData,isOutputData,dataType,numDims,dimArray,isFixedPoint,bias,slope,exponent,dataName,isComplex)\
 sf_debug_set_chart_data_props(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	(dataNumber),(dataScope),(isInputData),(isOutputData),\
	(dataType),(numDims),(dimArray),(isFixedPoint),(bias),(slope),(exponent),(dataName),(isComplex))
#define _SFD_STATE_INFO(v1,v2,v3)\
	sf_debug_set_chart_state_info(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1),(v2),(v3))
#define _SFD_CH_SUBSTATE_INDEX(v1,v2)\
	sf_debug_set_chart_substate_index(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1),(v2))
#define _SFD_ST_SUBSTATE_INDEX(v1,v2,v3)\
   sf_debug_set_chart_state_substate_index(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1),(v2),(v3))
#define _SFD_ST_SUBSTATE_COUNT(v1,v2)\
	sf_debug_set_chart_state_substate_count(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1),(v2))
#define _SFD_STATE_COV_WTS(v1,v2,v3,v4)\
	sf_debug_set_instance_state_coverage_weights(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,CHARTINSTANCE_INSTANCENUMBER,(v1),(v2),(v3),(v4))
#define _SFD_STATE_COV_MAPS(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10) \
 sf_debug_set_chart_state_coverage_maps(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
   (v1),(v2),(v3),(v4),(v5),(v6),(v7),(v8),(v9),(v10))
#define _SFD_TRANS_COV_WTS(v1,v2,v3,v4,v5) \
	sf_debug_set_instance_transition_coverage_weights(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,CHARTINSTANCE_INSTANCENUMBER,\
   (v1),(v2),(v3),(v4),(v5))
#define 	_SFD_TRANS_COV_MAPS(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13) \
	sf_debug_set_chart_transition_coverage_maps(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
      (v1),\
		(v2),(v3),(v4),\
		(v5),(v6),(v7),\
		(v8),(v9),(v10),\
		(v11),(v12),(v13))

#define _SFD_DATA_CHANGE_EVENT_COUNT(v1,v2) \
	sf_debug_set_number_of_data_with_change_event_for_chart(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	(v1),(v2))
#define _SFD_STATE_ENTRY_EVENT_COUNT(v1,v2) \
	sf_debug_set_number_of_states_with_entry_event_for_chart(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	(v1),(v2))
#define _SFD_STATE_EXIT_EVENT_COUNT(v1,v2) \
	sf_debug_set_number_of_states_with_exit_event_for_chart(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	(v1),(v2))
#define _SFD_EVENT_SCOPE(v1,v2)\
	sf_debug_set_chart_event_scope(_rov_design_analysisMachineNumber_,\
	CHARTINSTANCE_CHARTNUMBER,(v1),(v2))

#define _SFD_CH_SUBSTATE_COUNT(v1) \
	sf_debug_set_chart_substate_count(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1))
#define _SFD_CH_SUBSTATE_DECOMP(v1) \
	sf_debug_set_chart_decomposition(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,(v1))

#define _SFD_CV_INIT_CHART(v1,v2,v3,v4)\
 sf_debug_cv_init_chart(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	CHARTINSTANCE_INSTANCENUMBER,(v1),(v2),(v3),(v4))

#define _SFD_CV_INIT_STATE(v1,v2,v3,v4,v5,v6,v7,v8)\
	sf_debug_cv_init_state(_rov_design_analysisMachineNumber_,CHARTINSTANCE_CHARTNUMBER,\
	CHARTINSTANCE_INSTANCENUMBER,(v1),(v2),(v3),(v4),(v5),(v6),(v7),(v8))

#define _SFD_CV_INIT_TRANS(v1,v2,v3,v4,v5,v6)\
     sf_debug_cv_init_trans(_rov_design_analysisMachineNumber_,\
	  CHARTINSTANCE_CHARTNUMBER,\
	  CHARTINSTANCE_INSTANCENUMBER,\
	  (v1),(v2),(v3),(v4),(v5),(v6))
#endif

