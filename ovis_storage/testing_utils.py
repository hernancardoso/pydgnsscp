def count_rtcm_messages_by_id(messages):
    message_id_count = {}
    for message in messages:
        message_id = message['message_id']

        # Increment the count for this type_id
        if message_id in message_id_count:
            message_id_count[message_id] += 1
        else:
            message_id_count[message_id] = 1
    return message_id_count
