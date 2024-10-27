from pyrtcm.rtcmreader import RTCMReader

def decode_rtcm_message(message):
    decoded_data = {}
    parsed_rtcm = RTCMReader.parse(message, validate= 1)
    return parsed_rtcm.identity, parsed_rtcm